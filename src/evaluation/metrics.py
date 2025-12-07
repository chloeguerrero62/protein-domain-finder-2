"""
Evaluation metrics for domain detection
"""

import numpy as np
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def compute_domain_count_accuracy(n_pred, n_true):
    """Domain count accuracy metrics"""
    return {
        'exact_match': int(n_pred == n_true),
        'absolute_error': abs(n_pred - n_true),
        'relative_error': abs(n_pred - n_true) / n_true if n_true > 0 else 0
    }


def compute_boundary_f1(labels_pred, labels_true):
    """
    Compute F1 score for domain boundaries
    
    A boundary exists between residues i and i+1 if they have different labels
    """
    # Get boundaries
    pred_boundaries = set(i for i in range(len(labels_pred)-1) 
                         if labels_pred[i] != labels_pred[i+1])
    true_boundaries = set(i for i in range(len(labels_true)-1) 
                         if labels_true[i] != labels_true[i+1])
    
    if len(true_boundaries) == 0:
        return 0.0 if len(pred_boundaries) > 0 else 1.0
    
    # Compute precision, recall, F1
    tp = len(pred_boundaries & true_boundaries)
    fp = len(pred_boundaries - true_boundaries)
    fn = len(true_boundaries - pred_boundaries)
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    return f1


def compute_all_metrics(labels_pred, labels_true, n_true):
    """Compute all evaluation metrics"""
    n_pred = len(np.unique(labels_pred))
    
    metrics = {
        'n_predicted': n_pred,
        'n_true': n_true,
        **compute_domain_count_accuracy(n_pred, n_true)
    }
    
    if labels_true is not None:
        metrics.update({
            'ari': adjusted_rand_score(labels_true, labels_pred),
            'nmi': normalized_mutual_info_score(labels_true, labels_pred),
            'boundary_f1': compute_boundary_f1(labels_pred, labels_true)
        })
    
    return metrics
