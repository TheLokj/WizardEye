# -*- coding: utf-8 -*-

"""Project version metadata.

DISPLAY_VERSION is human-friendly and can be shown in CLI/help.
PACKAGE_VERSION is PEP 440 compliant for packaging tools.
"""

DISPLAY_VERSION = "beta-0.0.5"
PACKAGE_VERSION = "0.0.5b0"

import subprocess
import sys
from typing import Optional

def _get_git_commit_hash() -> Optional[str]:
	"""Get the current Git commit hash if available.
	
	Returns:
		Optional[str]: The Git commit hash as a string, or None if it cannot be determined.
	"""
	try:
		hash_bytes = subprocess.check_output(["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL)
		return hash_bytes.decode("utf-8").strip()
	except (subprocess.CalledProcessError, FileNotFoundError):
		return None
	
def print_version_message():
	"""Print a startup message with version and Git commit hash if available."""
	git_hash = _get_git_commit_hash()
	if git_hash and (sys.stdout.isatty() and sys.stderr.isatty()):
		print(f"{'-'*30} WizardEye v{PACKAGE_VERSION} \033[0;90m(#{git_hash[0:7]})\033[0m {'-'*30}")
	elif git_hash and (not sys.stdout.isatty() or not sys.stderr.isatty()):
		print(f"{'-'*30} WizardEye v{PACKAGE_VERSION} (#{git_hash[0:7]}) {'-'*30}")
	else:
		print(f"{'-'*30} WizardEye v{PACKAGE_VERSION} {'-'*30}\n")
