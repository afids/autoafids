import os

from appdirs import AppDirs


def get_download_dir():
    if "AUTOAFIDS_CACHE_DIR" in os.environ.keys():
        download_dir = os.environ["AUTOAFIDS_CACHE_DIR"]
    else:
        # create local download dir if it doesn't exist
        dirs = AppDirs("autoafids", "khanlab")
        download_dir = dirs.user_cache_dir
    return download_dir
