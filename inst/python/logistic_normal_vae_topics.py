#!/usr/bin/env python3
# VAE topic model for doc-term counts.
# Writes theta/phi per K and a model_metrics.csv with loglik/perplexity.

from __future__ import annotations

import argparse
import os
import sys
from typing import List

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import torch
    from torch import nn
    from torch.utils.data import Dataset, DataLoader
except Exception as exc:  # pragma: no cover
    sys.stderr.write("ERROR: torch is required for logistic-normal VAE training.\n")
    sys.stderr.write("Install with: pip install torch\n")
    sys.stderr.write(f"Details: {exc}\n")
    sys.exit(1)


class BowDataset(Dataset):
    def __init__(self, X_csr: sp.csr_matrix):
        self.X = X_csr

    def __len__(self) -> int:
        return self.X.shape[0]

    def __getitem__(self, idx: int) -> torch.Tensor:
        row = self.X.getrow(idx).toarray().ravel().astype(np.float32)
        return torch.from_numpy(row)


class LogisticNormalVAE(nn.Module):
    def __init__(self, vocab_size: int, n_topics: int, hidden: int = 128, dropout: float = 0.0):
        super().__init__()
        self.vocab_size = vocab_size
        self.n_topics = n_topics
        self.encoder = nn.Sequential(
            nn.Linear(vocab_size, hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        self.mu = nn.Linear(hidden, n_topics)
        self.logvar = nn.Linear(hidden, n_topics)
        self.beta = nn.Parameter(torch.randn(n_topics, vocab_size))

    def encode(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        h = self.encoder(torch.log1p(x))
        return self.mu(h), self.logvar(h)

    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        eps = torch.randn_like(mu)
        return mu + eps * torch.exp(0.5 * logvar)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        theta = torch.softmax(z, dim=1)
        phi = torch.softmax(self.beta, dim=1)
        return torch.matmul(theta, phi) + 1e-12

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z)
        return recon, mu, logvar

    def topic_word_dist(self) -> torch.Tensor:
        return torch.softmax(self.beta, dim=1)


class ETMLikeVAE(nn.Module):
    def __init__(self, vocab_size: int, n_topics: int, hidden: int = 128, dropout: float = 0.0, emb_dim: int | None = None):
        super().__init__()
        self.vocab_size = vocab_size
        self.n_topics = n_topics
        if emb_dim is None:
            emb_dim = hidden
        self.encoder = nn.Sequential(
            nn.Linear(vocab_size, hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        self.mu = nn.Linear(hidden, n_topics)
        self.logvar = nn.Linear(hidden, n_topics)
        self.topic_emb = nn.Parameter(torch.randn(n_topics, emb_dim))
        self.word_emb = nn.Parameter(torch.randn(vocab_size, emb_dim))

    def encode(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        h = self.encoder(torch.log1p(x))
        return self.mu(h), self.logvar(h)

    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        eps = torch.randn_like(mu)
        return mu + eps * torch.exp(0.5 * logvar)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        theta = torch.softmax(z, dim=1)
        beta = torch.matmul(self.topic_emb, self.word_emb.T)
        phi = torch.softmax(beta, dim=1)
        return torch.matmul(theta, phi) + 1e-12

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z)
        return recon, mu, logvar

    def topic_word_dist(self) -> torch.Tensor:
        beta = torch.matmul(self.topic_emb, self.word_emb.T)
        return torch.softmax(beta, dim=1)


class ModalitySplitVAE(nn.Module):
    def __init__(
        self,
        vocab_size: int,
        n_topics: int,
        gene_idx: np.ndarray,
        peak_idx: np.ndarray,
        hidden: int = 128,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.vocab_size = vocab_size
        self.n_topics = n_topics
        self.hidden = hidden
        self.register_buffer("gene_idx", torch.tensor(gene_idx, dtype=torch.long))
        self.register_buffer("peak_idx", torch.tensor(peak_idx, dtype=torch.long))
        self.gene_in = int(len(gene_idx))
        self.peak_in = int(len(peak_idx))
        self.gene_proj = nn.Linear(self.gene_in, hidden) if self.gene_in > 0 else None
        self.peak_proj = nn.Linear(self.peak_in, hidden) if self.peak_in > 0 else None
        self.fuse = nn.Sequential(
            nn.Linear(hidden * 2, hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        self.mu = nn.Linear(hidden, n_topics)
        self.logvar = nn.Linear(hidden, n_topics)
        self.beta = nn.Parameter(torch.randn(n_topics, vocab_size))

    def encode(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        batch_size = x.shape[0]
        if self.gene_in > 0:
            x_gene = torch.log1p(x.index_select(1, self.gene_idx))
            h_gene = torch.relu(self.gene_proj(x_gene))
        else:
            h_gene = torch.zeros((batch_size, self.hidden), device=x.device, dtype=x.dtype)
        if self.peak_in > 0:
            x_peak = torch.log1p(x.index_select(1, self.peak_idx))
            h_peak = torch.relu(self.peak_proj(x_peak))
        else:
            h_peak = torch.zeros((batch_size, self.hidden), device=x.device, dtype=x.dtype)
        h = torch.cat([h_gene, h_peak], dim=1)
        h = self.fuse(h)
        return self.mu(h), self.logvar(h)

    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        eps = torch.randn_like(mu)
        return mu + eps * torch.exp(0.5 * logvar)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        theta = torch.softmax(z, dim=1)
        phi = torch.softmax(self.beta, dim=1)
        return torch.matmul(theta, phi) + 1e-12

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z)
        return recon, mu, logvar

    def topic_word_dist(self) -> torch.Tensor:
        return torch.softmax(self.beta, dim=1)


def _parse_k_grid(text: str) -> List[int]:
    vals = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        vals.append(int(part))
    return sorted(set([v for v in vals if v > 1]))


def _prepare_matrix(doc_term: pd.DataFrame) -> tuple[sp.csr_matrix, List[str], List[str]]:
    doc_ids = doc_term["doc_id"].drop_duplicates().tolist()
    term_ids = doc_term["term_id"].drop_duplicates().tolist()
    doc_index = {d: i for i, d in enumerate(doc_ids)}
    term_index = {t: i for i, t in enumerate(term_ids)}
    rows = doc_term["doc_id"].map(doc_index).to_numpy()
    cols = doc_term["term_id"].map(term_index).to_numpy()
    data = doc_term["count"].astype(float).to_numpy()
    X = sp.csr_matrix((data, (rows, cols)), shape=(len(doc_ids), len(term_ids)))
    return X, doc_ids, term_ids


def _train_one(
    X: sp.csr_matrix,
    n_topics: int,
    epochs: int,
    batch_size: int,
    hidden: int,
    lr: float,
    seed: int,
    device: str,
    variant: str,
    term_ids: list[str],
    save_path: str | None = None,
) -> tuple[np.ndarray, np.ndarray, dict, str]:
    rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    if device == "cuda":
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
        else:
            device = "cpu"

    if variant == "vae_mlp":
        model = LogisticNormalVAE(X.shape[1], n_topics, hidden=hidden)
    elif variant == "moetm_encoder_decoder":
        model = ETMLikeVAE(X.shape[1], n_topics, hidden=hidden)
    elif variant == "multivi_encoder":
        terms = pd.Series(term_ids, dtype=str)
        gene_idx = np.where(terms.str.startswith("GENE:"))[0]
        peak_idx = np.where(terms.str.startswith("PEAK:"))[0]
        if gene_idx.size == 0 and peak_idx.size == 0:
            gene_idx = np.arange(len(term_ids))
        model = ModalitySplitVAE(X.shape[1], n_topics, gene_idx=gene_idx, peak_idx=peak_idx, hidden=hidden)
    else:
        raise ValueError(f"Unknown variant: {variant}")
    model = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    dataset = BowDataset(X)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, drop_last=False)

    model.train()
    for _epoch in range(int(epochs)):
        for batch in loader:
            batch = batch.to(device)
            recon, mu, logvar = model(batch)
            recon_loss = -(batch * torch.log(recon)).sum(dim=1)
            kl = 0.5 * torch.sum(torch.exp(logvar) + mu * mu - 1.0 - logvar, dim=1)
            loss = (recon_loss + kl).mean()
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

    model.eval()
    with torch.no_grad():
        # theta from encoder mean
        all_theta = []
        eval_loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, drop_last=False)
        for batch in eval_loader:
            batch = batch.to(device)
            mu, _logvar = model.encode(batch)
            theta = torch.softmax(mu, dim=1)
            all_theta.append(theta.cpu().numpy())
        theta = np.vstack(all_theta)
        phi = model.topic_word_dist().cpu().numpy()

        # metrics
        nll = 0.0
        n_tokens = float(X.sum())
        for batch in eval_loader:
            batch = batch.to(device)
            recon, mu, logvar = model(batch)
            recon_loss = -(batch * torch.log(recon)).sum(dim=1)
            nll += float(recon_loss.sum().cpu().numpy())

    if variant == "vae_mlp":
        variant_detail = "Encoder: shared MLP (log1p input)\nDecoder: free topic-term matrix (beta)"
    elif variant == "moetm_encoder_decoder":
        variant_detail = "Encoder: shared MLP (log1p input)\nDecoder: factorized topic_emb x word_emb (ETM-style)"
    elif variant == "multivi_encoder":
        variant_detail = "Encoder: modality-split MLP (gene/peak) with fusion\nDecoder: free topic-term matrix (beta)"
    else:
        variant_detail = "Encoder/Decoder: unspecified"

    model_summary = "\n".join(
        [
            "Model type: VAE topic model (neural network)",
            f"Variant: {variant}",
            f"Topics (K): {n_topics}",
            f"Vocab size: {X.shape[1]}",
            f"Hidden size: {hidden}",
            f"Device: {device}",
            "",
            variant_detail,
            "",
            "PyTorch module:",
            str(model),
        ]
    )

    if save_path:
        torch.save(
            {
                "state_dict": model.state_dict(),
                "config": {
                    "vocab_size": int(X.shape[1]),
                    "n_topics": int(n_topics),
                    "hidden": int(hidden),
                    "epochs": int(epochs),
                    "batch_size": int(batch_size),
                    "lr": float(lr),
                    "seed": int(seed),
                    "device": str(device),
                    "variant": str(variant),
                },
            },
            save_path,
        )

    perplexity = float(np.exp(nll / n_tokens)) if n_tokens > 0 else np.nan
    metrics = {
        "K": int(n_topics),
        "n_tokens": n_tokens,
        "nll": nll,
        "loglik": -nll,
        "perplexity": perplexity,
        "seed": int(seed),
        "epochs": int(epochs),
        "batch_size": int(batch_size),
        "hidden": int(hidden),
        "lr": float(lr),
        "variant": str(variant),
    }
    return theta, phi, metrics, model_summary


def main() -> None:
    ap = argparse.ArgumentParser(description="Logistic-normal VAE topic model for doc-term counts.")
    ap.add_argument("--doc-term", required=True, help="CSV with columns: doc_id, term_id, count")
    ap.add_argument("--out-dir", required=True, help="Output directory")
    ap.add_argument("--k-grid", required=True, help="Comma-separated K values, e.g. 5,10,15")
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--hidden", type=int, default=128)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--device", choices=["cpu", "cuda"], default="cpu")
    ap.add_argument(
        "--variant",
        choices=["vae_mlp", "moetm_encoder_decoder", "multivi_encoder"],
        default="vae_mlp",
    )
    args = ap.parse_args()

    doc_term = pd.read_csv(args.doc_term)
    if not {"doc_id", "term_id", "count"}.issubset(doc_term.columns):
        raise ValueError("doc_term must contain doc_id, term_id, count")
    doc_term = doc_term[doc_term["count"] > 0].copy()
    if doc_term.empty:
        raise ValueError("doc_term has no nonzero counts")

    X, doc_ids, term_ids = _prepare_matrix(doc_term)
    os.makedirs(args.out_dir, exist_ok=True)
    model_dir = os.path.join(args.out_dir, "vae_models")
    os.makedirs(model_dir, exist_ok=True)

    k_grid = _parse_k_grid(args.k_grid)
    if not k_grid:
        raise ValueError("K grid must include integers > 1")

    metrics_rows = []
    manifest_rows = []
    for k in k_grid:
        save_path = os.path.join(model_dir, f"model_K{k}.pt")
        theta, phi, metrics, model_summary = _train_one(
            X,
            n_topics=k,
            epochs=args.epochs,
            batch_size=args.batch_size,
            hidden=args.hidden,
            lr=args.lr,
            seed=args.seed,
            device=args.device,
            variant=args.variant,
            term_ids=term_ids,
            save_path=save_path,
        )
        metrics_rows.append(metrics)

        theta_df = pd.DataFrame(theta, index=doc_ids, columns=[f"Topic{i+1}" for i in range(k)])
        phi_df = pd.DataFrame(phi, index=[f"Topic{i+1}" for i in range(k)], columns=term_ids)
        theta_df.to_csv(os.path.join(model_dir, f"theta_K{k}.csv"))
        phi_df.to_csv(os.path.join(model_dir, f"phi_K{k}.csv"))
        arch_path = os.path.join(model_dir, f"model_K{k}_arch.txt")
        with open(arch_path, "w", encoding="utf-8") as handle:
            handle.write(model_summary)
            handle.write("\n")

        manifest_rows.append(
            {
                "variant": args.variant,
                "K": int(k),
                "model_path": save_path,
                "arch_path": arch_path,
                "theta_path": os.path.join(model_dir, f"theta_K{k}.csv"),
                "phi_path": os.path.join(model_dir, f"phi_K{k}.csv"),
                "vocab_size": int(X.shape[1]),
                "hidden": int(args.hidden),
                "epochs": int(args.epochs),
                "batch_size": int(args.batch_size),
                "lr": float(args.lr),
                "seed": int(args.seed),
                "device": str(args.device),
            }
        )

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(os.path.join(args.out_dir, "model_metrics.csv"), index=False)
    if manifest_rows:
        manifest_df = pd.DataFrame(manifest_rows)
        manifest_df.to_csv(os.path.join(args.out_dir, "vae_model_manifest.csv"), index=False)


if __name__ == "__main__":
    main()
