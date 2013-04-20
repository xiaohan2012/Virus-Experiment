def complex_pairwise_calc(complex_id_list, dataloader, calc_func, callback = None, exception_callback=None, symmetry=True):
    """
    (list of str, (str)=> compelx, (complex, complex) => double, (str, str, double) => ???) => ???

    performing complex pairwise calculation based on the calc_func it's provided, apply the callback if provided to the result
    also with the pairs that encounter exception during calculation
    
    """
    #get possible pairs
    from itertools import combinations, permutations, chain
    if symmetry:
        #a to b is the same as b to a
        pairs = chain(combinations(complex_id_list, 2), 
                      [(c,c) for c in complex_id_list])
    else:
        pairs = chain(permutations(complex_id_list, 2),
                      [(c,c) for c in complex_id_list])
        
    #calculate the pairwise value
    pairs_with_results = []
    exception_pairs = []
    for c1, c2 in pairs:
        try:
            pairs_with_results.append((c1, c2, calc_func(dataloader(c1),dataloader(c2))))
        except:
            exception_pairs.append((c1,c2))

    #if there is a callback, apply it to the results
    if callback:
        pairs_with_results = map(lambda tpls: callback(*tpls), pairs_with_results)
        
    if exception_callback:
        exception_pairs = map(lambda tpls: exception_callback(*tpls), exception_pairs)

    return pairs_with_results, exception_pairs
