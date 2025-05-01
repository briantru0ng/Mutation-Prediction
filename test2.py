import os
import ast
import json
from pathlib import Path

def extract_imports_from_code(code_str):
    try:
        tree = ast.parse(code_str)
        return {
            node.names[0].name if isinstance(node, ast.Import) else node.module.split('.')[0]
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom)) and node.module != '__future__'
        }
    except:
        return set()

def extract_from_py_and_ipynb(folder):
    imports = set()
    for file in Path(folder).rglob("*"):
        # print(file)
        if file.suffix == ".py":
            try:
                code = file.read_text(encoding='utf-8')
                imports |= extract_imports_from_code(code)
            except Exception as e:
                print(f"⚠️ Skipped {file}: {e}")
        elif file.suffix == ".ipynb":
            try:
                with open(file, encoding='utf-8') as f:
                    notebook = json.load(f)
                    for cell in notebook.get("cells", []):
                        if cell.get("cell_type") == "code":
                            cell_code = ''.join(cell.get("source", []))
                            imports |= extract_imports_from_code(cell_code)
            except Exception as e:
                print(f"⚠️ Skipped {file}: {e}")
    return sorted(imports)

# Generate requirements.txt
folder = ""
pkgs = extract_from_py_and_ipynb(folder)

with open("requirements.txt", "w") as f:
    f.write("\n".join(pkgs))

print("requirements.txt created with the following imports:")
print("\n".join(pkgs))
