import zipfile
import json
import tempfile
import shutil
import os

keras_path = "covid19_models/models/balanced_data_model.keras"
h5_path = "covid19_models/models/balanced_data_model.h5"
fixed_keras_path = "covid19_models/models/balanced_data_model_fixed.keras"

print(f"Fixing model file: {keras_path}")

def fix_batch_shape(obj):
    """Recursively fix batch_shape issues in config"""
    if isinstance(obj, dict):
        # Fix batch_shape -> shape
        if 'batch_shape' in obj and 'shape' not in obj:
            batch_shape = obj.pop('batch_shape')
            if batch_shape and len(batch_shape) > 1:
                obj['shape'] = batch_shape[1:]  # Remove batch dimension
            print(f"  Fixed batch_shape: {batch_shape} -> shape: {obj.get('shape')}")
        
        # Fix batch_input_shape -> input_shape
        if 'batch_input_shape' in obj and 'input_shape' not in obj:
            batch_input_shape = obj.pop('batch_input_shape')
            if batch_input_shape and len(batch_input_shape) > 1:
                obj['input_shape'] = batch_input_shape[1:]
            print(f"  Fixed batch_input_shape -> input_shape")
        
        # Recurse into nested dicts
        for key, value in list(obj.items()):
            fix_batch_shape(value)
    
    elif isinstance(obj, list):
        for item in obj:
            fix_batch_shape(item)

try:
    # Step 1: Extract and fix the model config
    with tempfile.TemporaryDirectory() as temp_dir:
        print("\nStep 1: Extracting model files...")
        
        # Extract all files
        with zipfile.ZipFile(keras_path, 'r') as zip_ref:
            zip_ref.extractall(temp_dir)
            print(f"  Extracted to: {temp_dir}")
            print(f"  Files: {zip_ref.namelist()}")
        
        # Read and fix config.json
        config_file = os.path.join(temp_dir, 'config.json')
        
        if os.path.exists(config_file):
            print("\nStep 2: Fixing config.json...")
            
            with open(config_file, 'r') as f:
                config = json.load(f)
            
            # Fix the config
            fix_batch_shape(config)
            
            # Write fixed config back
            with open(config_file, 'w') as f:
                json.dump(config, f, indent=2)
            
            print("  Config fixed!")
        else:
            print("  WARNING: config.json not found!")
        
        # Step 3: Create new zip with fixed config
        print("\nStep 3: Creating fixed model file...")
        
        with zipfile.ZipFile(fixed_keras_path, 'w', zipfile.ZIP_DEFLATED) as zip_out:
            for root, dirs, files in os.walk(temp_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    arcname = os.path.relpath(file_path, temp_dir)
                    zip_out.write(file_path, arcname)
        
        print(f"  Fixed model saved to: {fixed_keras_path}")
    
    # Step 4: Try to load the fixed model
    print("\nStep 4: Testing fixed model...")
    import tensorflow as tf
    
    try:
        model = tf.keras.models.load_model(fixed_keras_path, compile=False)
        print(f"✅ Fixed model loaded successfully!")
        print(f"   Layers: {len(model.layers)}")
        
        # Step 5: Save as .h5
        print("\nStep 5: Converting to .h5 format...")
        model.save(h5_path, save_format='h5')
        print(f"✅✅✅ Model saved to: {h5_path}")
        
        # Verify .h5 loads
        model_h5 = tf.keras.models.load_model(h5_path, compile=False)
        print(f"✅ .h5 model verified!")
        
        print("\n" + "="*70)
        print("SUCCESS! You can now use balanced_data_model.h5")
        print("="*70)
        
    except Exception as e:
        print(f"❌ Could not load fixed model: {e}")
        print("\nThe fixed .keras file is available at:")
        print(f"  {fixed_keras_path}")
        print("\nYou may need to try loading it in a different environment")

except Exception as e:
    print(f"❌ Error during fix process: {e}")
    import traceback
    traceback.print_exc()
