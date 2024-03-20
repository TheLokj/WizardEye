import math

def get_accuracy(TP, TN, FP, FN) :
    """
    This function calculates the accuracy thanks to the
    True Positive, True Negative, False Positive and False Negative data
    """
    return (TP+TN)/(TP+TN+FP+FN)

def get_f1score(TP, TN, FP, FN) :
    """
    This function calculates the F1-score thanks to the
    True Positive, True Negative, False Positive and False Negative data
    """
    return (2*TP)/(2*TP+FP+FN)

def get_mcc(TP, TN, FP, FN):
    """
    This function calculates the accuracy thanks to the
    True Positive, True Negative, False Positive and False Negative data
    """
    denom = math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    if denom == 0 :
        return None
    else :
        return (TP*TN - FP*FN)/denom