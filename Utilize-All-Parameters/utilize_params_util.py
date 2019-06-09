# Utility functions for use in all scripts here


def find_smirnoff_params() -> set:
    """Creates a set of all the smirnoff params"""
    smirnoff_ids = set()

    # number of parameters of each type
    num_params = {
        'b': 87,
        'a': 38,
        't': 158,
        'n': 35,
        'i': 4,
    }
    for (param_type, param_count) in num_params.items():
        for i in range(1, param_count + 1):
            smirnoff_ids.add(f"{param_type}{i}")

    return smirnoff_ids


def order_param_id(pid: str) -> (str, int):
    """Orders parameters by type then number"""
    return (pid[0], int(pid[1:]))
