def get_bool(value):
    val = value.strip().lower()
    true_values = {"true", "1", "yes", "on"}
    false_values = {"false", "0", "no", "off"}

    if val in true_values:
        return True
    elif val in false_values:
        return False
    else:
        return None