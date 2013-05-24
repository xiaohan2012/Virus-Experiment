def get_error_info(e):
    import traceback, os.path
    
    top = traceback.extract_stack()[-1]
    return ', '.join([type(e).__name__, os.path.basename(top[0]), str(top[1])])
