def calc_gc(sequence: str):
    gc_total = sum(1 for b in sequence if b in ('G', 'C', 'c', 'g'))
    return gc_total/len(sequence)