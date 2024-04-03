import logging


def debug(*args, **kwargs):
    if __debug__:
        logging.getLogger("feynamp").debug(*args, **kwargs)

def warning(*args, **kwargs):
    logging.getLogger("feynamp").warning(*args, **kwargs)

def info(*args, **kwargs):
    logging.getLogger("feynamp").info(*args, **kwargs)
