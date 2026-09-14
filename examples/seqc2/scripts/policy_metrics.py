"""Shared validation for development and held-out policy metrics."""

import math


def metric(value, label):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or not 0 <= value <= 1
    ):
        raise ValueError(f"{label} must be a finite number between 0 and 1")
    return float(value)


def f1(row, label):
    if all(key in row for key in ("tp", "fp", "fn")):
        tp, fp, fn = (float(row[key]) for key in ("tp", "fp", "fn"))
        if min(tp, fp, fn) < 0 or any(value != int(value) for value in (tp, fp, fn)):
            raise ValueError(f"{label} counts must be non-negative integers")
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        derived = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        if "f1" in row and abs(float(row["f1"]) - derived) > 1e-6:
            raise ValueError(f"{label} f1 disagrees with TP/FP/FN")
        return derived
    if "f1" in row:
        return metric(row["f1"], f"{label} f1")
    precision = metric(row.get("precision"), f"{label} precision")
    recall = metric(row.get("recall"), f"{label} recall")
    return 2 * precision * recall / (precision + recall) if precision + recall else 0.0
