from rcsbsearchapi.search import TextQuery
from rcsbsearchapi import rcsb_attributes as attrs

# OBTENCION DE PROTEINAS A PARTIR DE DATOS OBTENIDOS DE NCBI
gene_name = "SLC6A4"  # Nombre del gen a buscar
organism_name = "Homo sapiens"  # Nombre del organismo

# Consulta y obtencion de estructuras en PDB databank
q1 = attrs.rcsb_entity_source_organism.rcsb_gene_name.value == gene_name
q2 = attrs.rcsb_entity_source_organism.scientific_name == organism_name
query = q1 & q2
results = list(query())
print("Estructuras encontradas:", results)
#Si el nombre del gen no esta asociado a ninguna estructura dentro de PDB databank,
#el algortimo no generara resultados. 

#posterior analisis de las estructuras obtenidas, selecciona la que se trabajará con el resto del algoritmo
pdb_id = '5I6X'
print("Estructura seleccionada:",pdb_id)

import os  # Para manejar directorios
import requests  # Para descargar archivos

# Realiza un directorio para los archivos PDB 
protein_directory = "protein_structures"
#os.makedirs(protein_directory, exist_ok=True) #si no se tiene un directorio generado, desmarcar

# Descarga de la estructura PDB
pdb_request = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")

if pdb_request.status_code == 200:
    # Almacenamiento del archivo PDB
    with open(f"{protein_directory}/{pdb_id}.pdb", "w+") as f:
        f.write(pdb_request.text)
    print(f"Estructura {pdb_id} descargada y guardada en {protein_directory}/")
else:
    print(f"Error al descargar la estructura {pdb_id}. Código de estado: {pdb_request.status_code}")

#VISUALIZACION DE LA ESTRUCTURA DE PROTEINA
import MDAnalysis as mda
pdb_file = f"protein_structures/{pdb_id}.pdb"
u = mda.Universe(pdb_file)
print(u)
import nglview as nv
view=nv.show_mdanalysis(u) #visualizar estructura
view

#PREPARACION DEL BLANCO 
#Seleccion de atomos 
protein = u.select_atoms("protein")
ligand = u.select_atoms("resname 8PR") #el nombre del ligando va en funcion del ligando unido a la proteina seleccionada
water = u.select_atoms("resname HOH") #solo se llena en caso de que se quieran seleccionar las moleculas de agua
water 
 #Generacion de nuevo archivo .pdb
protein_pdb_path = f"{protein_directory}/protein_{pdb_id}.pdb"
protein.write(protein_pdb_path)
protein_pqr_path = f"{protein_directory}/protein_{pdb_id}.pqr"

# Comando para ejecutar PDB2PQR
import subprocess
pdb2pqr_command = [
    "pdb2pqr",
    "--ff=AMBER",  # Usar el campo de fuerza AMBER
    "--with-ph=7.0",  # Configurar pH
    protein_pdb_path,  # Archivo de entrada
    protein_pqr_path  # Archivo de salida
]

try:
    subprocess.run(pdb2pqr_command, check=True)
    print(f"Archivo PQR generado en: {protein_pqr_path}")
except subprocess.CalledProcessError as e:
    print(f"Error al ejecutar PDB2PQR: {e}")
    
pdbqt_directory = "pdbqt_structures" #directorio definido para los archivos pdbqt
os.makedirs(pdbqt_directory, exist_ok=True)  # Crea el directorio si no existe

# Carga la estructura PQR
pqr_path = f"{protein_directory}/protein_{pdb_id}.pqr"
u = mda.Universe(pqr_path)

pdbqt_path = f"{pdbqt_directory}/{pdb_id}.pdbqt" #conversion a PDBQT
u.atoms.write(pdbqt_path)

print(f"Archivo PDBQT generado y guardado en: {pdbqt_path}")

with open(f"{pdbqt_directory}/{pdb_id}.pdbqt", 'r') as file:
    file_content = file.read()
#MDAnalysis genera dos lineas que autodockVina no reconoce, por lo que es importante modificarlas 

# Reemplazar 'TITLE' y 'CRYST1' con 'REMARK'
file_content = file_content.replace('TITLE', 'REMARK').replace('CRYST1', 'REMARK')
with open(f"{pdbqt_directory}/{pdb_id}.pdbqt", 'w') as file:
    file.write(file_content)

#PREPARACION DE LOS LIGANDOS
#en este caso, los ligandos fueron previamente generados a traves de ChemDraw
#y minimizados molecularmente con MMimport os
import os
import subprocess
from pathlib import Path

# Define la ruta a la carpeta de entrada y salida
input_dir = Path(os.getcwd(), "Ligands")
output_dir = Path("C:/Users/User/Downloads/python/ligands_pdbqt")

# Verifica que la carpeta de entrada exista
if not input_dir.exists():
    print(f"Error: La carpeta {input_dir} no existe.")
    exit(1)

# Crea la carpeta de salida si no existe
output_dir.mkdir(parents=True, exist_ok=True)

# Lista todos los archivos en la carpeta de entrada
dir_list = os.listdir(input_dir)

# Itera sobre los archivos para convertir de .mol2 a .pdbqt
for molecule in dir_list:
    if molecule.endswith(".mol2"):
        # Rutas completas para el archivo de entrada y salida
        full_element_path = input_dir / molecule
        converted_name = molecule.replace(".mol2", ".pdbqt")
        full_converted_name = output_dir / converted_name

        # Ejecuta el comando de conversión usando Open Babel
        command = subprocess.run(
            ["obabel", str(full_element_path), "-O", str(full_converted_name)],
            text=True,
        )

        # Verifica si la conversión fue exitosa
        if command.returncode == 0:
            print(f"Convertido: {full_element_path} -> {full_converted_name}")
        else:
            print(f"Error al convertir {full_element_path}")
    else:
        print(f"Saltando archivo no compatible: {molecule}")

#DOCKING MOLECULAR
#encontrando coordenadas para la gridbox
import MDAnalysis as mda
#centro de gridbox
original_structure = mda.Universe("protein_structures/5I6X.pdb")
ligand_mda = original_structure.select_atoms("resname 8PR")
pocket_center =ligand_mda.center_of_geometry()
print(pocket_center)
#tamaño de gridbox
ligand_box = ligand_mda.positions.max(axis=0) -ligand_mda.positions.max(axis=0) +5  
ligand_box
#convirtiendo de Numpy arrays a listas
pocket_center = pocket_center.tolist()
ligand_box = ligand_box.tolist()
#docking con AutodockVina
import os
import csv
from pathlib import Path
from vina import Vina
pdb_id = "5I6X"
ligand = "8PR"
os.makedirs("Docking_results", exist_ok=True)
receptor_file = Path("C:/Users/User/Downloads/python/receptor.pdbqt")  # Archivo PDBQT del receptor
ligands_dir = Path("C:/Users/User/Downloads/python/ligands_pdbqt")  # Carpeta con ligandos
results_dir = Path("C:/Users/User/Downloads/python/docking_results")  # Carpeta para resultados de poses
csv_dir = Path("C:/Users/User/Downloads/python/docking_csv")  # Carpeta para archivos CSV
from vina import Vina
v = Vina(sf_name="vina")
v.set_receptor(f"pdbqt/{pdb_id}.pdbqt")
for ligand_file in ligands_dir.glob("*.pdbqt"):  # Esto solo seleccionará los archivos .pdbqt
    ligand_name = ligand_file.stem  # Obtener el nombre del ligando sin la extensión

    # Configurar el ligando
    v.set_ligand_from_file(str(ligand_file))  # Pasar la ruta completa del archivo de ligando
v.compute_vina_maps(center=pocket_center, box_size=ligand_box) 
print(f"Ejecutando docking para el ligando: {ligand_name}")
Vina.dock(exhaustiveness=8, n_poses=10)
# Guardar las poses en un archivo de salida
output_pose_file = results_dir / f"{ligand_name}_poses.pdbqt"
vina.write_poses(str(output_pose_file), n_poses=10)

# Obtener energías de acoplamiento
results = vina.poses()
docking_results = [
        {"receptor": pdb_id, "ligand": ligand_name, "pose": i + 1, "affinity": pose["score"]}
        for i, pose in enumerate(results)
    ]

 # Crear archivo CSV para el ligando en la carpeta docking_csv
output_csv_file = csv_dir / f"{ligand_name}_docking.csv"
with open(output_csv_file, mode="w", newline="") as csvfile:
        csv_writer = csv.DictWriter(csvfile, fieldnames=["receptor", "ligand", "pose", "affinity"])
        csv_writer.writeheader()
        csv_writer.writerows(docking_results)

print(f"Resultados guardados: {output_csv_file}")
