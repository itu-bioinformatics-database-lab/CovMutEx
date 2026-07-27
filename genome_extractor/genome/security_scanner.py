"""
Simple security scanner with file extension validation and ClamAV malware detection
"""

import os
from datetime import datetime

# Try to import ClamAV
try:
    import clamd
    CLAMAV_AVAILABLE = True
except ImportError:
    CLAMAV_AVAILABLE = False
    print("[WARNING] clamd not installed. Install with: pip install clamd")


# Allowed file extensions
ALLOWED_MODEL_EXTENSIONS = ['.keras', '.h5', '.hdf5', '.pkl', '.pt', '.pth', '.bin', '.onnx', '.safetensors']
ALLOWED_EXTRACTOR_EXTENSIONS = ['.py']
ALLOWED_HELPER_EXTENSIONS = ['.py', '.json', '.txt', '.csv', '.nwk', '.tsv', '.fasta', '.fa', '.yaml', '.yml', '.npy']


def validate_file_extension(filename, allowed_extensions):
    """
    Validate file extension against allowed list
    
    Args:
        filename: Name of the file
        allowed_extensions: List of allowed extensions (e.g., ['.py', '.json'])
    
    Returns:
        (is_valid: bool, message: str)
    """
    _, ext = os.path.splitext(filename.lower())
    
    if ext not in allowed_extensions:
        return False, f"Extension '{ext}' not allowed. Allowed: {', '.join(allowed_extensions)}"
    
    return True, f"Extension '{ext}' is valid"


def scan_file_for_malware(file_path):
    """
    Scan a single file for malware using ClamAV
    
    Args:
        file_path: Path to the file to scan
    
    Returns:
        (is_safe: bool, message: str)
    """
    if not CLAMAV_AVAILABLE:
        print(f"[WARNING] ClamAV not available - skipping malware scan for {file_path}")
        return True, "ClamAV not installed - scan skipped"
    
    try:
        # Try to connect to ClamAV daemon
        try:
            scanner = clamd.ClamdUnixSocket()
            scanner.ping()
        except:
            try:
                scanner = clamd.ClamdNetworkSocket()
                scanner.ping()
            except Exception as e:
                print(f"[WARNING] ClamAV daemon not running: {e}")
                return True, "ClamAV daemon not running - scan skipped"
        
        # Scan the file
        result = scanner.scan(file_path)
        
        if not result or file_path not in result:
            return True, "No threats detected"
        
        status, threat_name = result[file_path]
        
        if status == 'OK':
            return True, "No threats detected"
        else:
            return False, f"Malware detected: {threat_name}"
            
    except Exception as e:
        print(f"[ERROR] Malware scan failed: {e}")
        # Fail-safe: reject file if scan fails
        return False, f"Scan error: {str(e)}"


def scan_folder_for_malware(folder_path):
    """
    Scan all files in a folder for malware
    
    Args:
        folder_path: Path to folder containing uploaded files
    
    Returns:
        (all_clean: bool, results: dict)
    """
    results = {}
    all_clean = True
    
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        
        if os.path.isfile(file_path):
            is_safe, message = scan_file_for_malware(file_path)
            results[filename] = {
                "safe": is_safe,
                "message": message,
                "scanned_at": datetime.now().isoformat()
            }
            
            if not is_safe:
                all_clean = False
                print(f"[SECURITY] Malware detected in {filename}: {message}")
    
    return all_clean, results


def validate_uploaded_files(model_file=None, extractor_file=None, helper_files=None):
    """
    Validate file extensions for uploaded files
    
    Args:
        model_file: Model file object (optional)
        extractor_file: Feature extractor file object (optional)
        helper_files: Dict of helper file objects with custom names (optional)
    
    Returns:
        (is_valid: bool, error_message: str or None)
    """
    # Validate model file
    if model_file:
        is_valid, message = validate_file_extension(model_file.name, ALLOWED_MODEL_EXTENSIONS)
        if not is_valid:
            return False, f"Model file '{model_file.name}': {message}"
    
    # Validate extractor file
    if extractor_file:
        is_valid, message = validate_file_extension(extractor_file.name, ALLOWED_EXTRACTOR_EXTENSIONS)
        if not is_valid:
            return False, f"Extractor file '{extractor_file.name}': {message}"
    
    # Validate helper files
    if helper_files:
        # Handle dict, list, or set
        filenames = helper_files.values() if isinstance(helper_files, dict) else helper_files
        
        for custom_name in filenames:
            is_valid, message = validate_file_extension(custom_name, ALLOWED_HELPER_EXTENSIONS)
            if not is_valid:
                return False, f"Helper file '{custom_name}': {message}"
    
    return True, None
