def count_gc(seq: str) -> int:
    counter = 0
    for nucl in seq:
        if nucl.upper() == "G" or nucl.upper() == "C":
            counter += 1
    return counter


def is_filter_gc(seq: str, gc_bounds: tuple | float | int) -> bool:
    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)
    percentage_gc = count_gc(seq) / len(seq) * 100
    return gc_bounds[0] <= percentage_gc <= gc_bounds[1]


def is_filter_length(seq: str, length_bounds: tuple | int) -> bool:
    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)
    return length_bounds[0] <= len(seq) <= length_bounds[1]


def is_filter_quality(scores: str, quality_threshold: float | int) -> bool:
    q_scores = [(ord(score) - 33) for score in scores]
    avg_q_scores = sum(q_scores) / len(scores)
    return avg_q_scores >= quality_threshold
