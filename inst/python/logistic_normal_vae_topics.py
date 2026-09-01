#!/usr/bin/env python3
# VAE topic model for doc-term counts.
# Writes theta/phi per K and a model_metrics.csv with loglik/perplexity.

from __future__ import annotations

import argparse
import copy
import hashlib
import os
import sys
import time
import warnings
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


def _read_doc_term(path: str) -> pd.DataFrame:
    lower = path.lower()
    if lower.endswith((".arrow", ".feather", ".ipc")):
        try:
            return pd.read_feather(path)
        except Exception as exc:
            raise RuntimeError(f"Failed to read Arrow/Feather doc-term file: {path}; {exc}") from exc
    return pd.read_csv(path)


class BowDataset(Dataset):
    def __init__(self, X_csr: sp.csr_matrix):
        self.X = X_csr

    def __len__(self) -> int:
        return self.X.shape[0]

    def __getitem__(self, idx: int) -> int:
        return int(idx)

    def collate(self, indices: list[int]) -> torch.Tensor:
        rows = self.X[np.asarray(indices, dtype=np.int64), :]
        dense = rows.toarray().astype(np.float32, copy=False)
        return torch.from_numpy(dense)


def _bow_loader(
    dataset: BowDataset,
    batch_size: int,
    shuffle: bool,
) -> DataLoader:
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        drop_last=False,
        collate_fn=dataset.collate,
    )


class LogisticNormalVAE(nn.Module):
    def __init__(self, vocab_size: int, n_topics: int, hidden: int = 128, dropout: float = 0.0, topic_word_temperature: float = 1.0):
        super().__init__()
        self.vocab_size = vocab_size
        self.n_topics = n_topics
        self.topic_word_temperature = float(topic_word_temperature)
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
        logvar = torch.clamp(logvar, min=-10.0, max=10.0)
        eps = torch.randn_like(mu)
        return mu + eps * torch.exp(0.5 * logvar)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        theta = torch.softmax(z, dim=1)
        phi = torch.softmax(self.beta / self.topic_word_temperature, dim=1)
        return torch.matmul(theta, phi) + 1e-12

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z)
        return recon, mu, logvar

    def topic_word_dist(self) -> torch.Tensor:
        return torch.softmax(self.beta / self.topic_word_temperature, dim=1)


class ETMLikeVAE(nn.Module):
    def __init__(self, vocab_size: int, n_topics: int, hidden: int = 128, dropout: float = 0.0, emb_dim: int | None = None, topic_word_temperature: float = 1.0):
        super().__init__()
        self.vocab_size = vocab_size
        self.n_topics = n_topics
        self.topic_word_temperature = float(topic_word_temperature)
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
        logvar = torch.clamp(logvar, min=-10.0, max=10.0)
        eps = torch.randn_like(mu)
        return mu + eps * torch.exp(0.5 * logvar)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        theta = torch.softmax(z, dim=1)
        beta = torch.matmul(self.topic_emb, self.word_emb.T)
        phi = torch.softmax(beta / self.topic_word_temperature, dim=1)
        return torch.matmul(theta, phi) + 1e-12

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z)
        return recon, mu, logvar

    def topic_word_dist(self) -> torch.Tensor:
        beta = torch.matmul(self.topic_emb, self.word_emb.T)
        return torch.softmax(beta / self.topic_word_temperature, dim=1)


class ModalitySplitVAE(nn.Module):
    def __init__(
        self,
        vocab_size: int,
        n_topics: int,
        gene_idx: np.ndarray,
        peak_idx: np.ndarray,
        hidden: int = 128,
        dropout: float = 0.0,
        topic_word_temperature: float = 1.0,
    ):
        super().__init__()
        self.vocab_size = vocab_size
        self.n_topics = n_topics
        self.hidden = hidden
        self.topic_word_temperature = float(topic_word_temperature)
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
        logvar = torch.clamp(logvar, min=-10.0, max=10.0)
        eps = torch.randn_like(mu)
        return mu + eps * torch.exp(0.5 * logvar)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        theta = torch.softmax(z, dim=1)
        phi = torch.softmax(self.beta / self.topic_word_temperature, dim=1)
        return torch.matmul(theta, phi) + 1e-12

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z)
        return recon, mu, logvar

    def topic_word_dist(self) -> torch.Tensor:
        return torch.softmax(self.beta / self.topic_word_temperature, dim=1)


def _parse_k_grid(text: str) -> List[int]:
    vals = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        vals.append(int(part))
    return sorted(set([v for v in vals if v > 1]))


def _term_signature(term_ids: list[str]) -> str:
    digest = hashlib.sha256()
    for term_id in term_ids:
        digest.update(str(term_id).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def _append_progress(path: str | None, stage: str, k: int | str, status: str, message: str = "") -> None:
    if not path:
        return
    exists = os.path.exists(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a", encoding="utf-8") as handle:
        if not exists:
            handle.write("timestamp\tstage\tK\tstatus\tmessage\n")
        handle.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\t{stage}\t{k}\t{status}\t{message}\n")
        handle.flush()


def _cuda_available() -> bool:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return bool(torch.cuda.is_available())
        except Exception:
            return False


def _resolve_device(requested_device: str) -> str:
    if requested_device == "auto":
        return "cuda" if _cuda_available() else "cpu"
    if requested_device == "cuda" and not _cuda_available():
        return "cpu"
    return requested_device


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


def _largest_remainder_counts(score: np.ndarray, total: int) -> np.ndarray:
    score = np.asarray(score, dtype=float)
    if score.size == 0:
        return np.asarray([], dtype=np.float64)
    if total < score.size:
        raise ValueError("Model-token allocation cannot retain every nonzero term")
    score[~np.isfinite(score) | (score < 0)] = 0
    if not np.any(score > 0):
        score[:] = 1
    remaining = int(total - score.size)
    quota = score / score.sum() * remaining
    counts = np.ones(score.size, dtype=np.int64) + np.floor(quota).astype(np.int64)
    cumulative = np.cumsum(quota - np.floor(quota))
    rounded_cumulative = np.floor(cumulative + 1e-8).astype(np.int64)
    residual = np.diff(np.concatenate(([0], rounded_cumulative)))
    counts += residual
    correction = int(total - counts.sum())
    if correction:
        counts[int(np.argmax(counts))] += correction
    return counts.astype(np.float64)


def _finalize_condition_tf_counts(
    X: sp.csr_matrix,
    doc_ids: list[str],
    term_ids: list[str],
    peak_gene_ratio: float | None,
    condition_term_idf: bool,
    condition_term_idf_floor: float,
) -> tuple[sp.csr_matrix, dict]:
    if peak_gene_ratio is None and not condition_term_idf:
        return X, {
            "final_peak_gene_token_ratio": np.nan,
            "condition_term_idf": False,
            "condition_term_idf_floor": float(condition_term_idf_floor),
        }
    X = X.tocsr().astype(float)
    condition_ids = np.asarray([str(value).rsplit("::", 1)[0] for value in doc_ids])
    terms = np.asarray(term_ids, dtype=str)
    gene_term = np.char.startswith(terms, "GENE:")
    peak_term = np.char.startswith(terms, "PEAK:")
    if peak_gene_ratio is not None and np.any(~(gene_term | peak_term)):
        raise ValueError("Final Gene/Peak token balancing found unsupported terms")

    idf = np.ones(X.shape[1], dtype=float)
    unique_conditions = np.unique(condition_ids)
    if condition_term_idf and unique_conditions.size > 1:
        frequency = np.zeros(X.shape[1], dtype=np.int32)
        for condition in unique_conditions:
            rows = np.where(condition_ids == condition)[0]
            frequency += np.asarray(X[rows].getnnz(axis=0) > 0, dtype=np.int32)
        maximum_idf = np.log((unique_conditions.size + 1) / 2)
        specificity = np.log(
            (unique_conditions.size + 1) / (frequency + 1)
        ) / maximum_idf
        specificity = np.clip(specificity, 0, 1)
        idf = condition_term_idf_floor + (1 - condition_term_idf_floor) * specificity

    output_data = np.empty_like(X.data, dtype=float)
    source_tokens = 0.0
    final_tokens = 0.0
    peak_tokens = 0.0
    gene_tokens = 0.0
    for row in range(X.shape[0]):
        start, end = X.indptr[row], X.indptr[row + 1]
        columns = X.indices[start:end]
        values = X.data[start:end]
        scores = values * idf[columns]
        total = int(round(values.sum()))
        source_tokens += total
        if peak_gene_ratio is None:
            allocated = _largest_remainder_counts(scores, total)
        else:
            row_gene = gene_term[columns]
            row_peak = peak_term[columns]
            if not np.any(row_gene) or not np.any(row_peak):
                raise ValueError("Every condition::TF document must contain Gene and Peak terms")
            target_peak = int(round(total * peak_gene_ratio / (1 + peak_gene_ratio)))
            target_peak = max(int(row_peak.sum()), min(total - int(row_gene.sum()), target_peak))
            target_gene = total - target_peak
            allocated = np.empty(values.size, dtype=float)
            allocated[row_gene] = _largest_remainder_counts(scores[row_gene], target_gene)
            allocated[row_peak] = _largest_remainder_counts(scores[row_peak], target_peak)
        output_data[start:end] = allocated
        final_tokens += allocated.sum()
        gene_tokens += allocated[gene_term[columns]].sum()
        peak_tokens += allocated[peak_term[columns]].sum()
    output = sp.csr_matrix(
        (output_data, X.indices.copy(), X.indptr.copy()),
        shape=X.shape,
    )
    return output, {
        "final_peak_gene_token_ratio": (
            float(peak_gene_ratio) if peak_gene_ratio is not None else np.nan
        ),
        "condition_term_idf": bool(condition_term_idf),
        "condition_term_idf_floor": float(condition_term_idf_floor),
        "source_model_tokens": float(source_tokens),
        "final_model_tokens": float(final_tokens),
        "final_gene_tokens": float(gene_tokens),
        "final_peak_tokens": float(peak_tokens),
        "final_peak_gene_token_ratio_observed": (
            float(peak_tokens / gene_tokens) if gene_tokens > 0 else np.nan
        ),
    }


def _matched_gene_peak_indices(term_ids: list[str]) -> tuple[np.ndarray, np.ndarray]:
    term_position = {term: index for index, term in enumerate(term_ids)}
    genes = []
    peaks = []
    for term, gene_position in term_position.items():
        if not term.startswith("GENE:"):
            continue
        peak_position = term_position.get(f"PEAK:{term[5:]}")
        if peak_position is None:
            continue
        genes.append(gene_position)
        peaks.append(peak_position)
    return np.asarray(genes, dtype=np.int64), np.asarray(peaks, dtype=np.int64)


def _mapped_gene_peak_indices(
    term_ids: list[str],
    map_path: str | None,
) -> tuple[np.ndarray, np.ndarray]:
    exact_gene, exact_peak = _matched_gene_peak_indices(term_ids)
    if not map_path:
        return exact_gene, exact_peak
    mapping = pd.read_csv(map_path)
    if not {"peak_id", "gene_key"}.issubset(mapping.columns):
        raise ValueError("Peak-gene map must contain peak_id and gene_key")
    term_position = {term: index for index, term in enumerate(term_ids)}
    pair_set = set(zip(exact_gene.tolist(), exact_peak.tolist()))
    for row in mapping[["peak_id", "gene_key"]].itertuples(index=False):
        gene_position = term_position.get(f"GENE:{row.gene_key}")
        peak_position = term_position.get(f"PEAK:{row.peak_id}")
        if gene_position is not None and peak_position is not None:
            pair_set.add((gene_position, peak_position))
    if not pair_set:
        return np.asarray([], dtype=np.int64), np.asarray([], dtype=np.int64)
    ordered = sorted(pair_set)
    return (
        np.asarray([value[0] for value in ordered], dtype=np.int64),
        np.asarray([value[1] for value in ordered], dtype=np.int64),
    )


def _paired_term_weights(
    gene_idx: np.ndarray,
    weighting: str,
) -> np.ndarray:
    gene_idx = np.asarray(gene_idx, dtype=np.int64)
    if gene_idx.size == 0:
        return np.asarray([], dtype=np.float32)
    if weighting == "uniform":
        return np.ones(gene_idx.size, dtype=np.float32)
    if weighting != "target_equal":
        raise ValueError(f"Unknown paired-term weighting: {weighting}")
    _, inverse, counts = np.unique(
        gene_idx,
        return_inverse=True,
        return_counts=True,
    )
    return (1.0 / counts[inverse]).astype(np.float32)


def _split_validation_counts(
    X: sp.csr_matrix,
    fraction: float,
    seed: int,
) -> tuple[sp.csr_matrix, sp.csr_matrix | None]:
    if fraction <= 0:
        return X, None
    rng = np.random.default_rng(seed)
    validation = X.copy().astype(float)
    rounded = np.maximum(np.rint(validation.data).astype(np.int64), 0)
    held_out = rng.binomial(rounded, fraction)
    validation.data = held_out.astype(float)
    validation.eliminate_zeros()
    training = X.copy().astype(float)
    training.data = training.data - held_out
    training.data[training.data < 0] = 0
    training.eliminate_zeros()
    if training.nnz == 0 or validation.nnz == 0:
        return X, None
    return training.tocsr(), validation.tocsr()


def _validation_nll(
    model: nn.Module,
    training_dataset: Dataset,
    validation_dataset: Dataset,
    batch_size: int,
    device: str,
) -> float:
    value = 0.0
    model.eval()
    with torch.no_grad():
        for training_batch, validation_batch in zip(
            _bow_loader(training_dataset, batch_size, shuffle=False),
            _bow_loader(validation_dataset, batch_size, shuffle=False),
        ):
            training_batch = training_batch.to(device)
            validation_batch = validation_batch.to(device)
            mu, _ = model.encode(training_batch)
            recon = model.decode(mu)
            value += float(
                (-(validation_batch * torch.log(recon)).sum()).cpu().numpy()
            )
    return value


def _topic_hellinger_similarity(
    phi: torch.Tensor,
    term_idx: torch.Tensor,
) -> torch.Tensor:
    if phi.shape[0] < 2 or term_idx.numel() == 0:
        return phi.new_tensor(0.0)
    modality_phi = phi.index_select(1, term_idx)
    modality_phi = modality_phi / modality_phi.sum(dim=1, keepdim=True).clamp_min(1e-12)
    rooted = torch.sqrt(modality_phi.clamp_min(1e-12))
    similarity = torch.matmul(rooted, rooted.T)
    off_diagonal = similarity.sum() - torch.diagonal(similarity).sum()
    return off_diagonal / float(phi.shape[0] * (phi.shape[0] - 1))


def _topic_diversity_loss(
    phi: torch.Tensor,
    gene_idx: torch.Tensor,
    peak_idx: torch.Tensor,
) -> torch.Tensor:
    components = []
    if gene_idx.numel() > 0:
        components.append(_topic_hellinger_similarity(phi, gene_idx))
    if peak_idx.numel() > 0:
        components.append(_topic_hellinger_similarity(phi, peak_idx))
    if not components:
        return phi.new_tensor(0.0)
    return torch.stack(components).mean()


def _document_topic_separation_loss(mu: torch.Tensor) -> torch.Tensor:
    if mu.shape[1] < 2:
        return mu.new_tensor(0.0)
    theta = torch.softmax(mu, dim=1).clamp_min(1e-12)
    conditional_entropy = -(theta * theta.log()).sum(dim=1).mean()
    marginal = theta.mean(dim=0).clamp_min(1e-12)
    marginal_entropy = -(marginal * marginal.log()).sum()
    return (conditional_entropy - marginal_entropy) / np.log(theta.shape[1])


def _document_anchor_topic_probabilities(
    X: sp.csr_matrix,
    n_topics: int,
    seed: int,
) -> np.ndarray:
    counts = X.astype(np.float64, copy=False).tocsr()
    totals = np.maximum(np.asarray(counts.sum(axis=1)).ravel(), 1.0)
    profiles = counts.multiply(1.0 / totals[:, None]).tocsr()
    rooted = profiles.copy()
    rooted.data = np.sqrt(rooted.data)
    center = np.asarray(rooted.mean(axis=0)).ravel()
    row_norm = np.asarray(rooted.multiply(rooted).sum(axis=1)).ravel()
    first = int(np.argmax(row_norm - 2.0 * np.asarray(rooted @ center).ravel()))
    anchors = [first]
    while len(anchors) < min(n_topics, profiles.shape[0]):
        selected = rooted[np.asarray(anchors)]
        similarities = (rooted @ selected.T).toarray()
        min_distance = np.min(1.0 - similarities, axis=1)
        min_distance[np.asarray(anchors)] = -np.inf
        anchors.append(int(np.argmax(min_distance)))

    corpus = np.asarray(counts.sum(axis=0)).ravel()
    corpus = corpus / np.maximum(corpus.sum(), 1.0)
    anchor_profiles = profiles[np.asarray(anchors)].toarray()
    topic_profiles = [
        0.99 * anchor_profiles[index] + 0.01 * corpus
        for index in range(anchor_profiles.shape[0])
    ]
    if n_topics > len(topic_profiles):
        rng = np.random.default_rng(seed)
        for _ in range(n_topics - len(topic_profiles)):
            weights = rng.dirichlet(np.full(profiles.shape[0], 0.3))
            mixture = np.asarray(weights @ profiles).ravel()
            topic_profiles.append(0.99 * mixture + 0.01 * corpus)
    out = np.vstack(topic_profiles)
    out = np.maximum(out, 1e-12)
    return out / out.sum(axis=1, keepdims=True)


def _initialize_topics_from_documents(
    model: nn.Module,
    X: sp.csr_matrix,
    n_topics: int,
    seed: int,
) -> bool:
    if not hasattr(model, "beta"):
        return False
    probabilities = _document_anchor_topic_probabilities(X, n_topics, seed)
    logits = np.log(probabilities)
    logits = logits - logits.mean(axis=1, keepdims=True)
    with torch.no_grad():
        model.beta.copy_(
            torch.tensor(logits, dtype=model.beta.dtype, device=model.beta.device)
        )
    return True


def _load_warm_start_model(
    model: nn.Module,
    path: str,
    *,
    vocab_size: int,
    n_topics: int,
    hidden: int,
    variant: str,
    topic_word_temperature: float,
    term_ids: list[str],
    allow_temperature_change: bool = False,
) -> int:
    if not os.path.isfile(path):
        raise ValueError(f"Warm-start model does not exist: {path}")
    checkpoint = torch.load(path, map_location="cpu", weights_only=True)
    if not isinstance(checkpoint, dict) or "state_dict" not in checkpoint:
        raise ValueError("Warm-start model must contain state_dict and config")
    config = checkpoint.get("config", {})
    expected = {
        "vocab_size": int(vocab_size),
        "n_topics": int(n_topics),
        "hidden": int(hidden),
        "variant": str(variant),
    }
    if not allow_temperature_change:
        expected["topic_word_temperature"] = float(topic_word_temperature)
    mismatches = [
        f"{name}={config.get(name)!r} (expected {value!r})"
        for name, value in expected.items()
        if config.get(name) != value
    ]
    if mismatches:
        raise ValueError(
            "Warm-start model is incompatible: " + "; ".join(mismatches)
        )
    checkpoint_term_signature = config.get("term_signature")
    if (
        checkpoint_term_signature is not None
        and checkpoint_term_signature != _term_signature(term_ids)
    ):
        raise ValueError("Warm-start model uses a different ordered vocabulary")
    try:
        model.load_state_dict(checkpoint["state_dict"], strict=True)
    except RuntimeError as exc:
        raise ValueError(f"Warm-start model tensor shapes are incompatible: {exc}") from exc
    return int(config.get("total_epochs_completed", config.get("epochs_completed", 0)))


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
    paired_term_regularization: float = 0.0,
    paired_term_weighting: str = "uniform",
    topic_diversity_regularization: float = 0.0,
    document_topic_separation_regularization: float = 0.0,
    topic_initialization: str = "random",
    topic_word_temperature: float = 1.0,
    peak_gene_map: str | None = None,
    validation_fraction: float = 0.0,
    validation_interval: int = 1,
    early_stopping_patience: int = 0,
    early_stopping_relative_min_delta: float = 1e-3,
    lr_reduction_patience: int = 0,
    lr_reduction_factor: float = 0.5,
    minimum_lr: float = 1e-6,
    weight_decay: float = 0.0,
    kl_warmup_epochs: int = 0,
    dropout: float = 0.0,
    warm_start_model: str | None = None,
    allow_warm_start_temperature_change: bool = False,
    save_path: str | None = None,
) -> tuple[np.ndarray, np.ndarray, dict, str]:
    rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    requested_device = str(device)
    device = _resolve_device(requested_device)
    if device == "cuda":
        torch.cuda.manual_seed_all(seed)

    if variant == "vae_mlp":
        model = LogisticNormalVAE(
            X.shape[1], n_topics, hidden=hidden, dropout=dropout,
            topic_word_temperature=topic_word_temperature,
        )
    elif variant == "moetm_encoder_decoder":
        model = ETMLikeVAE(
            X.shape[1], n_topics, hidden=hidden, dropout=dropout,
            topic_word_temperature=topic_word_temperature,
        )
    elif variant == "multivi_encoder":
        terms = pd.Series(term_ids, dtype=str)
        gene_idx = np.where(terms.str.startswith("GENE:"))[0]
        peak_idx = np.where(terms.str.startswith("PEAK:"))[0]
        if gene_idx.size == 0 and peak_idx.size == 0:
            gene_idx = np.arange(len(term_ids))
        model = ModalitySplitVAE(
            X.shape[1],
            n_topics,
            gene_idx=gene_idx,
            peak_idx=peak_idx,
            hidden=hidden,
            dropout=dropout,
            topic_word_temperature=topic_word_temperature,
        )
    else:
        raise ValueError(f"Unknown variant: {variant}")
    model = model.to(device)
    warm_start_epochs = 0
    topic_initialization_applied = False
    if warm_start_model:
        warm_start_epochs = _load_warm_start_model(
            model,
            warm_start_model,
            vocab_size=X.shape[1],
            n_topics=n_topics,
            hidden=hidden,
            variant=variant,
            topic_word_temperature=topic_word_temperature,
            term_ids=term_ids,
            allow_temperature_change=allow_warm_start_temperature_change,
        )
    elif topic_initialization == "document_anchor":
        topic_initialization_applied = _initialize_topics_from_documents(
            model,
            X,
            n_topics,
            seed,
        )
    optimizer = torch.optim.Adam(
        model.parameters(), lr=lr, weight_decay=weight_decay
    )
    scheduler = None
    if validation_fraction > 0 and lr_reduction_patience > 0:
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer,
            mode="min",
            factor=lr_reduction_factor,
            patience=lr_reduction_patience,
            min_lr=minimum_lr,
        )
    X_train, X_validation = _split_validation_counts(
        X, validation_fraction, seed + 104729
    )
    dataset = BowDataset(X_train)
    loader = _bow_loader(dataset, batch_size, shuffle=True)
    pair_gene_idx, pair_peak_idx = _mapped_gene_peak_indices(
        term_ids, peak_gene_map
    )
    pair_weight = _paired_term_weights(pair_gene_idx, paired_term_weighting)
    pair_gene_idx_t = torch.tensor(pair_gene_idx, dtype=torch.long, device=device)
    pair_peak_idx_t = torch.tensor(pair_peak_idx, dtype=torch.long, device=device)
    pair_weight_t = torch.tensor(pair_weight, dtype=torch.float32, device=device)
    terms = pd.Series(term_ids, dtype=str)
    diversity_gene_idx_t = torch.tensor(
        np.where(terms.str.startswith("GENE:"))[0],
        dtype=torch.long,
        device=device,
    )
    diversity_peak_idx_t = torch.tensor(
        np.where(terms.str.startswith("PEAK:"))[0],
        dtype=torch.long,
        device=device,
    )

    validation_dataset = (
        BowDataset(X_validation) if X_validation is not None else None
    )
    best_state = None
    best_validation_nll = np.inf
    epochs_without_improvement = 0
    epochs_completed = 0
    early_stopping_start_epoch = max(1, int(kl_warmup_epochs))
    if warm_start_model and validation_dataset is not None:
        best_validation_nll = _validation_nll(
            model,
            dataset,
            validation_dataset,
            batch_size,
            device,
        )
        best_state = copy.deepcopy(model.state_dict())
        print(
            "training "
            f"K={n_topics} epoch={warm_start_epochs} "
            f"validation_nll={best_validation_nll:.6g} "
            "warm_start_baseline=true",
            flush=True,
        )
    model.train()
    for _epoch in range(int(epochs)):
        epochs_completed = _epoch + 1
        for batch in loader:
            batch = batch.to(device)
            recon, mu, logvar = model(batch)
            logvar = torch.clamp(logvar, min=-10.0, max=10.0)
            recon_loss = -(batch * torch.log(recon)).sum(dim=1)
            kl = 0.5 * torch.sum(torch.exp(logvar) + mu * mu - 1.0 - logvar, dim=1)
            kl_scale = 1.0
            if kl_warmup_epochs > 0:
                kl_scale = min(1.0, (_epoch + 1) / float(kl_warmup_epochs))
            loss = (recon_loss + kl_scale * kl).mean()
            if paired_term_regularization > 0 and pair_gene_idx_t.numel() > 0:
                phi = model.topic_word_dist()
                gene_topic = phi.index_select(1, pair_gene_idx_t)
                peak_topic = phi.index_select(1, pair_peak_idx_t)
                gene_topic = gene_topic / gene_topic.sum(dim=0, keepdim=True).clamp_min(1e-12)
                peak_topic = peak_topic / peak_topic.sum(dim=0, keepdim=True).clamp_min(1e-12)
                mean_topic = 0.5 * (gene_topic + peak_topic)
                pair_js_by_pair = 0.5 * (
                    (
                        gene_topic
                        * (
                            gene_topic.clamp_min(1e-12).log()
                            - mean_topic.clamp_min(1e-12).log()
                        )
                    ).sum(dim=0)
                    + (
                        peak_topic
                        * (
                            peak_topic.clamp_min(1e-12).log()
                            - mean_topic.clamp_min(1e-12).log()
                        )
                    ).sum(dim=0)
                )
                pair_js = (
                    pair_js_by_pair * pair_weight_t
                ).sum() / pair_weight_t.sum().clamp_min(1e-12)
                token_scale = batch.sum(dim=1).mean().detach().clamp_min(1.0)
                loss = (
                    loss
                    + float(paired_term_regularization) * token_scale * pair_js
                )
            if topic_diversity_regularization > 0 and n_topics > 1:
                phi = model.topic_word_dist()
                diversity_loss = _topic_diversity_loss(
                    phi,
                    diversity_gene_idx_t,
                    diversity_peak_idx_t,
                )
                token_scale = batch.sum(dim=1).mean().detach().clamp_min(1.0)
                loss = (
                    loss
                    + float(topic_diversity_regularization)
                    * token_scale
                    * diversity_loss
                )
            if document_topic_separation_regularization > 0 and n_topics > 1:
                separation_loss = _document_topic_separation_loss(mu)
                token_scale = batch.sum(dim=1).mean().detach().clamp_min(1.0)
                loss = (
                    loss
                    + float(document_topic_separation_regularization)
                    * token_scale
                    * separation_loss
                )
            if not torch.isfinite(loss):
                raise RuntimeError(
                    f"Non-finite VAE loss for variant={variant}, K={n_topics}, epoch={_epoch + 1}"
                )
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()

        validate_now = (
            validation_dataset is not None
            and (
                epochs_completed % validation_interval == 0
                or epochs_completed == int(epochs)
            )
        )
        if validate_now:
            validation_nll = _validation_nll(
                model,
                dataset,
                validation_dataset,
                batch_size,
                device,
            )
            if scheduler is not None:
                scheduler.step(validation_nll)
            early_stopping_ready = epochs_completed >= early_stopping_start_epoch
            if early_stopping_ready:
                if (
                    epochs_completed == early_stopping_start_epoch
                    and not warm_start_model
                ):
                    best_validation_nll = np.inf
                    best_state = None
                    epochs_without_improvement = 0
                required_improvement = max(
                    1e-6,
                    abs(best_validation_nll)
                    * float(early_stopping_relative_min_delta),
                )
                if (
                    not np.isfinite(best_validation_nll)
                    or best_validation_nll - validation_nll
                    > required_improvement
                ):
                    best_validation_nll = validation_nll
                    best_state = copy.deepcopy(model.state_dict())
                    epochs_without_improvement = 0
                else:
                    epochs_without_improvement += 1
                if (
                    epochs_completed == early_stopping_start_epoch
                    or epochs_completed % 10 == 0
                ):
                    print(
                        "training "
                        f"K={n_topics} epoch={warm_start_epochs + epochs_completed} "
                        f"validation_nll={validation_nll:.6g} "
                        f"best_validation_nll={best_validation_nll:.6g} "
                        f"epochs_without_improvement={epochs_without_improvement}",
                        flush=True,
                    )
            model.train()
            if (
                early_stopping_ready
                and early_stopping_patience > 0
                and epochs_without_improvement >= early_stopping_patience
            ):
                break

    if best_state is not None:
        model.load_state_dict(best_state)

    model.eval()
    with torch.no_grad():
        # theta from encoder mean
        all_theta = []
        inference_dataset = BowDataset(X)
        eval_loader = _bow_loader(inference_dataset, batch_size, shuffle=False)
        for batch in eval_loader:
            batch = batch.to(device)
            mu, _logvar = model.encode(batch)
            theta = torch.softmax(mu, dim=1)
            all_theta.append(theta.cpu().numpy())
        theta = np.vstack(all_theta)
        phi = model.topic_word_dist().cpu().numpy()
        if not np.isfinite(theta).all() or not np.isfinite(phi).all():
            raise RuntimeError(f"Non-finite VAE theta/phi for variant={variant}, K={n_topics}")
        if not np.any(theta > 0) or not np.any(phi > 0):
            raise RuntimeError(f"Empty VAE theta/phi for variant={variant}, K={n_topics}")

        # metrics
        nll = 0.0
        n_tokens = float(X.sum())
        for batch in eval_loader:
            batch = batch.to(device)
            recon, mu, logvar = model(batch)
            recon_loss = -(batch * torch.log(recon)).sum(dim=1)
            nll += float(recon_loss.sum().cpu().numpy())

    paired_term_js = np.nan
    paired_term_argmax_agreement = np.nan
    if pair_gene_idx.size > 0:
        gene_topic = phi[:, pair_gene_idx]
        peak_topic = phi[:, pair_peak_idx]
        gene_topic = gene_topic / np.maximum(gene_topic.sum(axis=0, keepdims=True), 1e-12)
        peak_topic = peak_topic / np.maximum(peak_topic.sum(axis=0, keepdims=True), 1e-12)
        mean_topic = 0.5 * (gene_topic + peak_topic)
        pair_js_by_pair = 0.5 * (
            np.sum(
                gene_topic
                * (
                    np.log(np.maximum(gene_topic, 1e-12))
                    - np.log(np.maximum(mean_topic, 1e-12))
                ),
                axis=0,
            )
            + np.sum(
                peak_topic
                * (
                    np.log(np.maximum(peak_topic, 1e-12))
                    - np.log(np.maximum(mean_topic, 1e-12))
                ),
                axis=0,
            )
        )
        paired_term_js = float(
            np.sum(pair_js_by_pair * pair_weight)
            / np.maximum(np.sum(pair_weight), 1e-12)
        )
        pair_agreement = (
            np.argmax(gene_topic, axis=0)
            == np.argmax(peak_topic, axis=0)
        ).astype(float)
        paired_term_argmax_agreement = float(
            np.sum(pair_agreement * pair_weight)
            / np.maximum(np.sum(pair_weight), 1e-12)
        )

    phi_tensor = torch.tensor(phi, dtype=torch.float32)
    topic_hellinger_gene = float(
        _topic_hellinger_similarity(
            phi_tensor,
            torch.tensor(
                np.where(pd.Series(term_ids, dtype=str).str.startswith("GENE:"))[0],
                dtype=torch.long,
            ),
        ).item()
    )
    topic_hellinger_peak = float(
        _topic_hellinger_similarity(
            phi_tensor,
            torch.tensor(
                np.where(pd.Series(term_ids, dtype=str).str.startswith("PEAK:"))[0],
                dtype=torch.long,
            ),
        ).item()
    )
    theta_safe = np.maximum(theta, 1e-12)
    theta_conditional_entropy = float(
        np.mean(-np.sum(theta_safe * np.log(theta_safe), axis=1))
        / np.log(theta.shape[1])
    )
    theta_marginal = np.maximum(theta.mean(axis=0), 1e-12)
    theta_marginal_entropy = float(
        -np.sum(theta_marginal * np.log(theta_marginal))
        / np.log(theta.shape[1])
    )
    theta_mutual_information = float(
        theta_marginal_entropy - theta_conditional_entropy
    )

    if variant == "vae_mlp":
        variant_detail = "Encoder: shared MLP (log1p input)\nDecoder: free topic-term matrix (beta)"
    elif variant == "moetm_encoder_decoder":
        variant_detail = "Encoder: shared MLP (log1p input)\nDecoder: factorized topic_emb x word_emb (ETM-style)"
    elif variant == "multivi_encoder":
        variant_detail = "\n".join(
            [
                "Encoder: modality-split MLP (gene/peak) with fusion",
                "Decoder: free topic-term matrix (beta)",
                (
                    "Matched Gene/Peak topic regularization: "
                    f"{paired_term_regularization}"
                ),
            ]
        )
    else:
        variant_detail = "Encoder/Decoder: unspecified"

    model_summary = "\n".join(
        [
            "Model type: VAE topic model (neural network)",
            f"Variant: {variant}",
            f"Topics (K): {n_topics}",
            f"Vocab size: {X.shape[1]}",
            f"Hidden size: {hidden}",
            f"Requested device: {requested_device}",
            f"Resolved device: {device}",
            f"Warm-start model: {warm_start_model or 'none'}",
            f"Warm-start epochs: {warm_start_epochs}",
            f"Total epochs completed: {warm_start_epochs + epochs_completed}",
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
                    "term_signature": _term_signature(term_ids),
                    "n_topics": int(n_topics),
                    "hidden": int(hidden),
                    "epochs": int(epochs),
                    "epochs_completed": int(epochs_completed),
                    "warm_start_model": str(warm_start_model or ""),
                    "warm_start_epochs": int(warm_start_epochs),
                    "total_epochs_completed": int(warm_start_epochs + epochs_completed),
                    "batch_size": int(batch_size),
                    "lr": float(lr),
                    "seed": int(seed),
                    "requested_device": str(requested_device),
                    "resolved_device": str(device),
                    "variant": str(variant),
                    "paired_term_regularization": float(paired_term_regularization),
                    "paired_term_weighting": str(paired_term_weighting),
                    "topic_diversity_regularization": float(topic_diversity_regularization),
                    "document_topic_separation_regularization": float(document_topic_separation_regularization),
                    "topic_initialization": str(topic_initialization),
                    "topic_initialization_applied": bool(topic_initialization_applied),
                    "topic_word_temperature": float(topic_word_temperature),
                    "matched_gene_peak_terms": int(pair_gene_idx.size),
                    "validation_fraction": float(validation_fraction),
                    "validation_interval": int(validation_interval),
                    "early_stopping_patience": int(early_stopping_patience),
                    "early_stopping_relative_min_delta": float(
                        early_stopping_relative_min_delta
                    ),
                    "lr_reduction_patience": int(lr_reduction_patience),
                    "lr_reduction_factor": float(lr_reduction_factor),
                    "minimum_lr": float(minimum_lr),
                    "weight_decay": float(weight_decay),
                    "kl_warmup_epochs": int(kl_warmup_epochs),
                    "dropout": float(dropout),
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
        "epochs_completed": int(epochs_completed),
        "warm_start_model": str(warm_start_model or ""),
        "warm_start_epochs": int(warm_start_epochs),
        "total_epochs_completed": int(warm_start_epochs + epochs_completed),
        "batch_size": int(batch_size),
        "hidden": int(hidden),
        "lr": float(lr),
        "variant": str(variant),
        "requested_device": str(requested_device),
        "resolved_device": str(device),
        "paired_term_regularization": float(paired_term_regularization),
        "paired_term_weighting": str(paired_term_weighting),
        "topic_diversity_regularization": float(topic_diversity_regularization),
        "document_topic_separation_regularization": float(document_topic_separation_regularization),
        "topic_initialization": str(topic_initialization),
        "topic_initialization_applied": bool(topic_initialization_applied),
        "topic_word_temperature": float(topic_word_temperature),
        "topic_hellinger_gene": topic_hellinger_gene,
        "topic_hellinger_peak": topic_hellinger_peak,
        "theta_conditional_entropy": theta_conditional_entropy,
        "theta_marginal_entropy": theta_marginal_entropy,
        "theta_mutual_information": theta_mutual_information,
        "matched_gene_peak_terms": int(pair_gene_idx.size),
        "paired_term_js": paired_term_js,
        "paired_term_argmax_agreement": paired_term_argmax_agreement,
        "validation_fraction": float(validation_fraction),
        "validation_interval": int(validation_interval),
        "validation_nll": (
            float(best_validation_nll)
            if np.isfinite(best_validation_nll)
            else np.nan
        ),
        "early_stopping_patience": int(early_stopping_patience),
        "early_stopping_relative_min_delta": float(
            early_stopping_relative_min_delta
        ),
        "lr_reduction_patience": int(lr_reduction_patience),
        "lr_reduction_factor": float(lr_reduction_factor),
        "minimum_lr": float(minimum_lr),
        "final_lr": float(optimizer.param_groups[0]["lr"]),
        "weight_decay": float(weight_decay),
        "kl_warmup_epochs": int(kl_warmup_epochs),
        "dropout": float(dropout),
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
    ap.add_argument("--paired-term-regularization", type=float, default=0.0)
    ap.add_argument(
        "--paired-term-weighting",
        choices=["uniform", "target_equal"],
        default="uniform",
    )
    ap.add_argument("--topic-diversity-regularization", type=float, default=0.0)
    ap.add_argument("--document-topic-separation-regularization", type=float, default=0.0)
    ap.add_argument(
        "--topic-initialization",
        choices=["random", "document_anchor"],
        default="random",
    )
    ap.add_argument("--topic-word-temperature", type=float, default=1.0)
    ap.add_argument("--peak-gene-map", default=None)
    ap.add_argument("--validation-fraction", type=float, default=0.0)
    ap.add_argument("--validation-interval", type=int, default=1)
    ap.add_argument("--early-stopping-patience", type=int, default=0)
    ap.add_argument(
        "--early-stopping-relative-min-delta", type=float, default=1e-3
    )
    ap.add_argument("--lr-reduction-patience", type=int, default=0)
    ap.add_argument("--lr-reduction-factor", type=float, default=0.5)
    ap.add_argument("--minimum-lr", type=float, default=1e-6)
    ap.add_argument("--weight-decay", type=float, default=0.0)
    ap.add_argument("--kl-warmup-epochs", type=int, default=0)
    ap.add_argument("--dropout", type=float, default=0.0)
    ap.add_argument("--final-peak-gene-token-ratio", type=float, default=None)
    ap.add_argument("--condition-term-idf", action="store_true")
    ap.add_argument("--condition-term-idf-floor", type=float, default=0.1)
    ap.add_argument("--count-divisor", type=float, default=1.0)
    ap.add_argument(
        "--warm-start-model",
        default=None,
        help="Optional compatible model_K*.pt whose network weights initialize training",
    )
    ap.add_argument(
        "--allow-warm-start-temperature-change",
        action="store_true",
        help="Allow an intentional temperature change when tensor shapes match",
    )
    ap.add_argument("--device", choices=["cpu", "cuda", "auto"], default="auto")
    ap.add_argument(
        "--variant",
        choices=["vae_mlp", "moetm_encoder_decoder", "multivi_encoder"],
        default="vae_mlp",
    )
    ap.add_argument("--progress-log", default=None, help="Optional TSV progress log path")
    args = ap.parse_args()
    if not np.isfinite(args.paired_term_regularization) or args.paired_term_regularization < 0:
        raise ValueError("--paired-term-regularization must be a finite non-negative number")
    if not np.isfinite(args.topic_diversity_regularization) or args.topic_diversity_regularization < 0:
        raise ValueError("--topic-diversity-regularization must be a finite non-negative number")
    if not np.isfinite(args.document_topic_separation_regularization) or args.document_topic_separation_regularization < 0:
        raise ValueError("--document-topic-separation-regularization must be a finite non-negative number")
    if not np.isfinite(args.topic_word_temperature) or args.topic_word_temperature <= 0:
        raise ValueError("--topic-word-temperature must be finite and positive")
    if not np.isfinite(args.validation_fraction) or not 0 <= args.validation_fraction < 1:
        raise ValueError("--validation-fraction must be in [0, 1)")
    if args.validation_interval < 1:
        raise ValueError("--validation-interval must be >= 1")
    if args.early_stopping_patience < 0:
        raise ValueError("--early-stopping-patience must be non-negative")
    if (
        not np.isfinite(args.early_stopping_relative_min_delta)
        or args.early_stopping_relative_min_delta < 0
    ):
        raise ValueError(
            "--early-stopping-relative-min-delta must be non-negative"
        )
    if not np.isfinite(args.weight_decay) or args.weight_decay < 0:
        raise ValueError("--weight-decay must be non-negative")
    if args.lr_reduction_patience < 0:
        raise ValueError("--lr-reduction-patience must be non-negative")
    if not np.isfinite(args.lr_reduction_factor) or not 0 < args.lr_reduction_factor < 1:
        raise ValueError("--lr-reduction-factor must be in (0, 1)")
    if not np.isfinite(args.minimum_lr) or args.minimum_lr <= 0:
        raise ValueError("--minimum-lr must be positive")
    if args.kl_warmup_epochs < 0:
        raise ValueError("--kl-warmup-epochs must be non-negative")
    if not np.isfinite(args.dropout) or not 0 <= args.dropout < 1:
        raise ValueError("--dropout must be in [0, 1)")
    if (
        args.final_peak_gene_token_ratio is not None
        and (
            not np.isfinite(args.final_peak_gene_token_ratio)
            or args.final_peak_gene_token_ratio <= 0
        )
    ):
        raise ValueError("--final-peak-gene-token-ratio must be positive")
    if (
        not np.isfinite(args.condition_term_idf_floor)
        or not 0 <= args.condition_term_idf_floor <= 1
    ):
        raise ValueError("--condition-term-idf-floor must be in [0, 1]")
    if not np.isfinite(args.count_divisor) or args.count_divisor < 1:
        raise ValueError("--count-divisor must be finite and >= 1")

    try:
        torch_threads = int(os.environ.get("OMP_NUM_THREADS", os.environ.get("NUMEXPR_NUM_THREADS", "1")))
    except ValueError:
        torch_threads = 1
    if torch_threads > 0:
        torch.set_num_threads(torch_threads)

    resolved_device = _resolve_device(args.device)
    _append_progress(
        args.progress_log,
        "device",
        "all",
        "resolved",
        f"requested_device={args.device};resolved_device={resolved_device};cuda_available={_cuda_available()}",
    )

    _append_progress(args.progress_log, "read_doc_term", "all", "start", args.doc_term)
    doc_term = _read_doc_term(args.doc_term)
    if not {"doc_id", "term_id", "count"}.issubset(doc_term.columns):
        raise ValueError("doc_term must contain doc_id, term_id, count")
    doc_term = doc_term[doc_term["count"] > 0].copy()
    if doc_term.empty:
        raise ValueError("doc_term has no nonzero counts")

    X, doc_ids, term_ids = _prepare_matrix(doc_term)
    X, input_weighting = _finalize_condition_tf_counts(
        X,
        doc_ids,
        term_ids,
        args.final_peak_gene_token_ratio,
        args.condition_term_idf,
        args.condition_term_idf_floor,
    )
    if args.count_divisor > 1:
        X.data = np.maximum(1.0, np.rint(X.data / args.count_divisor))
        X.eliminate_zeros()
    input_weighting["count_divisor"] = float(args.count_divisor)
    _append_progress(
        args.progress_log,
        "prepare_matrix",
        "all",
        "done",
        f"docs={X.shape[0]};terms={X.shape[1]};nnz={X.nnz}",
    )
    os.makedirs(args.out_dir, exist_ok=True)
    model_dir = os.path.join(args.out_dir, "vae_models")
    os.makedirs(model_dir, exist_ok=True)

    k_grid = _parse_k_grid(args.k_grid)
    if not k_grid:
        raise ValueError("K grid must include integers > 1")

    metrics_rows = []
    manifest_rows = []
    for k in k_grid:
        t_k = time.time()
        _append_progress(args.progress_log, "train_k", int(k), "start", f"variant={args.variant}")
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
            paired_term_regularization=args.paired_term_regularization,
            paired_term_weighting=args.paired_term_weighting,
            topic_diversity_regularization=args.topic_diversity_regularization,
            document_topic_separation_regularization=args.document_topic_separation_regularization,
            topic_initialization=args.topic_initialization,
            topic_word_temperature=args.topic_word_temperature,
            peak_gene_map=args.peak_gene_map,
            validation_fraction=args.validation_fraction,
            validation_interval=args.validation_interval,
            early_stopping_patience=args.early_stopping_patience,
            early_stopping_relative_min_delta=(
                args.early_stopping_relative_min_delta
            ),
            lr_reduction_patience=args.lr_reduction_patience,
            lr_reduction_factor=args.lr_reduction_factor,
            minimum_lr=args.minimum_lr,
            weight_decay=args.weight_decay,
            kl_warmup_epochs=args.kl_warmup_epochs,
            dropout=args.dropout,
            warm_start_model=args.warm_start_model,
            allow_warm_start_temperature_change=(
                args.allow_warm_start_temperature_change
            ),
            save_path=save_path,
        )
        metrics_rows.append(metrics)
        metrics.update(input_weighting)
        _append_progress(
            args.progress_log,
            "train_k",
            int(k),
            "model_done",
            f"elapsed_sec={round(time.time() - t_k, 1)};loglik={metrics.get('loglik')}",
        )

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
                "requested_device": str(args.device),
                "resolved_device": str(metrics.get("resolved_device")),
                "paired_term_regularization": float(args.paired_term_regularization),
                "validation_fraction": float(args.validation_fraction),
                "validation_nll": metrics.get("validation_nll"),
                "epochs_completed": int(metrics.get("epochs_completed", args.epochs)),
                "warm_start_model": str(args.warm_start_model or ""),
                "warm_start_epochs": int(metrics.get("warm_start_epochs", 0)),
                "total_epochs_completed": int(
                    metrics.get("total_epochs_completed", args.epochs)
                ),
                "early_stopping_relative_min_delta": float(
                    args.early_stopping_relative_min_delta
                ),
                "weight_decay": float(args.weight_decay),
                "kl_warmup_epochs": int(args.kl_warmup_epochs),
                "dropout": float(args.dropout),
                "matched_gene_peak_terms": int(metrics.get("matched_gene_peak_terms", 0)),
                "paired_term_js": metrics.get("paired_term_js"),
                "paired_term_argmax_agreement": metrics.get(
                    "paired_term_argmax_agreement"
                ),
            }
        )
        _append_progress(args.progress_log, "train_k", int(k), "done", f"elapsed_sec={round(time.time() - t_k, 1)}")

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(os.path.join(args.out_dir, "model_metrics.csv"), index=False)
    if manifest_rows:
        manifest_df = pd.DataFrame(manifest_rows)
        manifest_df.to_csv(os.path.join(args.out_dir, "vae_model_manifest.csv"), index=False)
    _append_progress(args.progress_log, "all", "all", "done", f"n_models={len(metrics_rows)}")


if __name__ == "__main__":
    main()
