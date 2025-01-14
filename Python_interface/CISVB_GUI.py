#####################################################################################################################################################
# python GUI interface of Chemical Intuitive & Symmetric Valence Bond structure set generation code to interact with the users;
# This interface can be used to get control keywords (Class CtrlInput): Chemical Formula (need to validate other inputs & create name of output file) 
# , Number of Active orbitals, Number of Active Electrons, Multiplicity; Geometry of molecule (Class ReadGeo) : there are two options to insert
# geometry (1: can be brows to upload a file, 2: can insert manually), The geometry can be inserted in two different units Bohr or Angstrom;
# Orbital informations (class Orb_Input): User need to specify the number of the active atoms (according to geometry) associated with active orbitals.
# User also need to specify type of the orbitals (sig or pi(x) or pi(y) or pi(z)), they can also provide fragments (the cluster of atomic orbitals 
# such as s, px, py, pz, dxx, dxy, dyy etc. over which the orbitals delocalised); finally user need to specify spatial keywords (class Keywd_Input)
# : method type (Chem_Inst: provide Chemically meaningfull VB structure set or Rumer: provide Rumer VB structure set), Rumer Set type (Single_Rumer_set:
# Rumer set corresponding to the present orbital order or All_Rumer_Set: Rumer sets corresponding to all possible unique permutations of the orbital 
# order); Structure Type (Covalent or Ionic or Both); Chem Inst Set Type (Single set: One Best Chemically meaning ful set, All Best Sets: All highesr 
# quality sets, All Sets: all possible sets with their quality scores, Eq Bond: Sets where all bonds appear same number of time ); Chem_Inst Str Type 
# ( Symmetry: sets will be symmetric, Asymmetry: there may not be any symmetry in the sets, Checksymm: Check if the set or sets are symmetric or not);
# Symmetric group arrangement (Quality-Arrange: Symmetric groups are arranged according to chemical quality scores, Big-to-Small: big group to small 
# group, small-to-big: arrange small group to big group); Estimate Overlap (Yes: Estimate the overlap of the structures of a set or No:); Decide 
# priorities of the Chemical Qualities: (The Chemical quality of a structures has been measured depending of 5 criteria which are, 1) IAB : Presence of
# inter atomic bonds in the structure -Negative quality, as less as this type of bonds prefered 2) NAB: Near atomic Bonds - positive quality -as many 
# as this type of bond is prefered, 3) SBB: Symmetry Breaking Bonds - negative quality, PDB: Pre-defined bonds, user can provide their prefered bonds,
# taken as positive quality, PDR: Pre-defined radicals, for systems with radicals, user can provide there prefered radical positions or orbitals, taken
# as positive quality), User can choose the priorities of these qualities; Priority Structures (user can put some structures they prefered to have in 
# the set, it is possible only if the provided structures are linearly independent); Maximum Number of Output File: if All sets are selected then it 
# might be important for larger system to specify output set numbers, by default one output set can have only 75000 sets.
# After getting all the data, class Run_Fort send them to Fortran subroutines to generate set of VB structures and store them in the output file.      
########################################################################################################################################################
    
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
#import tkinter.messagebox as msg
import re
import os
#import program_main
import subprocess
import sys
#import cisvb
import numpy as np
import Pmw
import symm_str

class GlobVar:
    "Class Contains Global Variables"
    num_orbital= None           # number of active localised VB orbitals taken in the calculation.
    num_electron = None         # number of active lectrons taken in the calculations.
    multiplicity = None         # multiplicity = 2*S + 1/2 , S = Total Spin of the system.
    type_orb_count = None       # orbital can be sigma, px, py or pz; type_orb_count counts how many different type of orbs present.
    num_iao = None              # number of inactive orbitals in the system.
    geometry_inserted = False   # This indicates that the geometry has not been provided.
    readgeo = None              # Instance of ReadGeo class, it indicates if the instance has been created or not.
    at_list_bold = [                                                                                       
        'H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE', 'NA', 'MG', 'AL',                   # List of Atoms according to periodic table
        'SL', 'P', 'S', 'CL', 'AR', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NR',
        'CU', 'ZN', 'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC',
        'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', 'CE',
        'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA',
        'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA'
    ]


class Ctrl_Input:
    # Class create the input panes to get Controll keywords: Chemical Formula (need to validate other inputs & create name of output file) 
    # Number of Active orbitals, Number of Active Electrons, Multiplicity
    def __init__(self, root):
        self.root = root
        self.input_text=''
        self.insert = False
        self.orbital_button = None
        self.unit_type_entry = False
        self.file_path = None

        Pmw.initialise(root)
        self.balloon = Pmw.Balloon(root)
        style_colour_frame = ttk.Style()
        style_colour_frame.configure("Colour_Frame.TFrame", background="lightblue")
        style_colour_label = ttk.Style()
        style_colour_label.configure("Colour_Label.TLabel", foreground="black", background="Ivory", relief="flat", font=("Arial", 14))
        style_colour_label1 = ttk.Style()
        style_colour_label1.configure("Colour_Label1.TLabel", foreground="black", background="lightblue", font=("Arial", 14))
        style_colour_label2 = ttk.Style()
        style_colour_label2.configure("Colour_Label2.TLabel", foreground="black", background="deepskyblue", font=("Arial", 14))
        style = ttk.Style()
        style.configure("Custom.TRadiobutton", foreground="black", background="Lightblue", relief="raised", font=("Arial", 14))

        # Main Frame
        self.frame = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.frame1 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
        self.frame1.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.frame2 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
        self.frame2.grid(row=2, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.frame3 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
        self.frame3.grid(row=3, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.frame4 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
        self.frame4.grid(row=4, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.rect_border = ttk.Frame(self.frame2, relief="solid", borderwidth=2, padding=5)
        self.rect_border.grid(row=0, column=1)

        self.unit_type = tk.StringVar(value="None")
        # read prefered file name
        self.read_filename()
        self.create_geo_unit()

        # Input Fields
        self.entries = {}
        self.keywords = [" Number of Active Orbitals (nao) : ", " Number of Active Electrons (nae) : ", " Multiplicity of Active Part (nmul) : "]
        self.create_ctrl_pans()

        # read geometry from a file
        brows_geo_lebel = ttk.Label(self.frame1, text = "Brows to Upload Geometry File", style = "Colour_Label.TLabel")
        brows_geo_lebel.grid(row = 0, column=0, sticky=tk.W, padx=5, pady=5)
        self.balloon.bind(brows_geo_lebel, "If you have the geometry saved in any dat file you can brows that file \n" 
                          "and insert the geometry here using 'Brows' button; in the geometry file you \n" 
                          "should have 4 or 5 columns first column: Atoms, second column: atomic numbers,\n" 
                          "third column: x coordinates, fourth column: y coordinates, fifth column: z coordinates" )

        brows_geo_button = ttk.Button(self.frame1, text = "Brows", command = self.read_geometry)
        brows_geo_button.grid(row = 0, column = 1, sticky = tk.W, padx = 5, pady = 5)
        self.balloon.bind(brows_geo_button, "If you have the geometry saved in any dat file you can brows that file \n" 
                          "and insert the geometry here using 'Brows' button; in the geometry file you \n" 
                          "should have 4 or 5 columns first column: Atoms, second column: atomic numbers,\n" 
                          "third column: x coordinates, fourth column: y coordinates, fifth column: z coordinates" )

        # upload geometry manually
        manual_geo_lebel = ttk.Label(self.frame1, text = "Insert Geometry Manually", style = "Colour_Label.TLabel")
        manual_geo_lebel.grid(row = 2, column=0, sticky=tk.W, padx=5, pady=5)
        self.balloon.bind(manual_geo_lebel, "If you wish to insert the geometry manually please click 'Geometry' button \n" 
                          " in the geometry file you should have 4 or 5 columns first column: Atoms, second column: atomic numbers,\n" 
                          "third column: x coordinates, fourth column: y coordinates, fifth column: z coordinates" )

        manual_geo_button = ttk.Button(self.frame1, text = "Geometry", command = self.insert_geo_manually)
        manual_geo_button.grid(row = 2, column = 1, sticky = tk.W, padx = 5, pady = 5)
        self.balloon.bind(manual_geo_button, "If you wish to insert the geometry manually please click 'Geometry' button \n" 
                          " in the geometry file you should have 4 or 5 columns first column: Atoms, second column: atomic numbers,\n" 
                          "third column: x coordinates, fourth column: y coordinates, fifth column: z coordinates" )

        insert_button = ttk.Button(self.frame4, text="Insert", command=self.validate_and_generate)
        insert_button.grid(row=4, column=1, padx=5, pady=10)
        self.balloon.bind(insert_button, "Press Insert after providing all control keywords: Without inserting \n"
                          "these control inputs you are not able to go to the next step")


    def create_geo_unit(self):
        # Its creats the radio button to get the information about the units of geometry provided.
        # there are only two options: Bohr or Angstrom. 
        units = ["Bohr", "Angs"]
        label = ttk.Label(self.frame1, text = "Unit of the Geometry Data", style = "Colour_Label.TLabel")
        label.grid(row = 3, column = 0, sticky = tk.W, padx = 10, pady = 10)
        for i, unit in enumerate(units, start=1):
            button = ttk.Radiobutton(
                self.frame1,
                text=unit,
                value=unit,
                variable=self.unit_type,
                command=self.update_geo_unit,
                style="Custom.TRadiobutton"
            )
            button.grid(row=3, column=i, padx=10, pady=10)

    def update_geo_unit(self):
        geo_unit = self.unit_type.get()
        if geo_unit:
            self.unit_type_entry = True
            print('geo_unit',geo_unit)
            return (geo_unit)

    def read_filename(self):
        molecule = {}
        label1 = ttk.Label(self.frame, text="Please enter the chemical formula", style = "Colour_Label.TLabel")
        label1.grid(row = 0, column=0, sticky=tk.W, padx=5, pady=5)
        label2 = ttk.Label(self.frame, text="(e.g.: C6H6)", style = "Colour_Label1.TLabel")
        label2.grid(row = 0, column=3, sticky=tk.W, padx=5, pady=5)
        self.balloon.bind(label1, "Enter a sequence of element symbols followed by numbers \n" 
                          "to specify the amounts of desired elements (e.g., C6H6, N2O,) \n"
                          "For reaction use 'under_score' in between two reactants (e.g., O_H2)") 
        self.molecule_entry = ttk.Entry(self.frame, width=15)
        self.molecule_entry.grid(row = 0, column=1, padx=5, pady=5)

    def count_Total_Electron(self, mol_entry):
        molecule = {}
        try:
            pattern = re.compile(r"([A-Z][a-z]?)(\d*)")
            # Remove underscores for handling O_H2 as OH2
            input_string = mol_entry.replace("_", "")
            # Find all matches in the string
            matches = pattern.findall(input_string)
            for element, count in matches:
                # If no count is given, assume it's 1
                count = int(count) if count else 1
                # Add element to dictionary or update count if already present
                molecule[element] = molecule.get(element, 0) + count
            total = 0
            for element, count in molecule.items():
                if element in GlobVar.at_list_bold:
                    atomic_number = GlobVar.at_list_bold.index(element) + 1  # Atomic number = index + 1
                    total += atomic_number * count
                else:
                    raise ValueError(f"Element {element} is not spelled Correctly.")
            return total
        except ValueError:
            messagebox.showerror("Incorrect Molecular Formula","The molecular formula is not in Correct Format \n"
                                     "Enter a sequence of element symbols followed by numbers \n" 
                          "to specify the amounts of desired elements (e.g., C6H6, N2O,) \n"
                          "For reaction use 'under_score' in between two reactants (e.g., O_H2)") 
            return



    def read_geometry(self):
        # Open file dialog to select a file
#        global readgeo
        self.file_path = filedialog.askopenfilename(title="Select Geometry File")

        if self.file_path:  # Check if a file is selected
            GlobVar.readgeo = Read_Geo(self.file_path)  # Initialize the Read_Geo class
            GlobVar.readgeo.read_geometry()       # Call the read_geometry method

            ttk.Button(
                self.frame1,
                text=f"View_Geometry",
                command=GlobVar.readgeo.display_geometry
            ).grid(row=0, column=2, columnspan=2)



    def insert_geo_manually(self):
        """Allows manual insertion of geometry data via Read_Geo."""
#        global readgeo
        GlobVar.readgeo = Read_Geo(self.file_path)  # Initialize the Read_Geo class
        GlobVar.readgeo.insert_geo(self.root)  # Call insert_geo method
        if GlobVar.geometry_inserted == True:
            button = ttk.Button(self.frame1, text = "View Geometry", command = GlobVar.readgeo.display_geometry)
            button.grid(row = 2, column = 2, sticky = tk.W, padx = 5, pady = 5)



    def create_ctrl_pans(self):
        lebel2 = ttk.Label(self.rect_border, text="Enter Control (ctrl) Keywords", font=("Arial", 14))
        lebel2.grid(row = 2, column=0, columnspan=2, pady=5)

        for idx, key in enumerate(self.keywords):
            ttk.Label(self.frame4, text=key, style = "Colour_Label.TLabel").grid(row = idx + 1, column=0, sticky=tk.W, padx=5, pady=5)
            entry = ttk.Entry(self.frame4, width=20)
            entry.grid(row = idx + 1, column=1, padx=5, pady=5)
            self.entries[key] = entry

    def validate_and_generate(self):
#        global num_iao, geometry_inserted, num_orbital, num_electron, multiplicity
        self.ctrl_inputs = []
        if GlobVar.geometry_inserted == False:
            messagebox.showerror("Geometry Error","Dont forget to insert the geometry of the system")
            
        try:
            molecule_string = self.molecule_entry.get()
            if not molecule_string:
                messagebox.showerror("Molecule Error","Please enter molecular Formula")
                return
            end_name = "_structures.dat"
            self.output_file_name =f"{molecule_string}{end_name}"

            Total_Electrons = self.count_Total_Electron(molecule_string)
            print('Total Number of Electrons', Total_Electrons)
        except ValueError:
            messagebox.showerror("Molecule Error","Please enter molecular Formula Properly ")

        if self.unit_type_entry == False:
            messagebox.showerror("Unit Type Missing", "Please Select a Unit between Bohr and Angs (angstrom).")
            return
        else:
            # Loop through each input field to validate
            for key, entry in self.entries.items():
                value = entry.get().strip()
                if not value:  # If any field is empty
                    messagebox.showerror("Validation Error", f"'{key}' is empty. Please fill in all fields.")
                    entry.configure(background="pink")  # Highlight the empty field
                    entry.focus_set()  # Focus the empty field
                    self.insert = False
                    return False
                else:
                    self.ctrl_inputs.append(value)

            for entry in self.entries.values():
                entry.configure(background="white")  # Reset background color
            self.insert = True
            if self.orbital_button:  # Enable the Orbital button
                nao, nae, nmul = self.ctrl_inputs
                GlobVar.num_orbital= int(nao)
                GlobVar.num_electron= int(nae)
                GlobVar.multiplicity= int(nmul)
                print('nao, nae, nmul',GlobVar.num_orbital, GlobVar.num_electron, GlobVar.multiplicity)
                if (Total_Electrons- GlobVar.num_electron) % 2== 0:
                    GlobVar.num_iao = int((Total_Electrons- GlobVar.num_electron)/2)
                    print('num_iao',GlobVar.num_iao)
                else:
                    messagebox.showerror("Active Orbital Error","The number of Active Orbitals or \n"
                                         "Number of Active Electrons or Molecular Formula is not correct")
                    return
                self.orbital_button.config(state=tk.NORMAL)
            return True  # Validation successful

#    def generate_ctrl_input(self):
#        if not self.insert:
#            messagebox.showerror("Validation Error", "Please insert valid values before insert other data.")
#            return
#       # return(nao, nae, nmul)


    def get_ctrl_keywds(self):
        geometry_unit = self.update_geo_unit()
        nao, nae, nmul= self.ctrl_inputs
        return(geometry_unit, nao, nae, nmul, self.output_file_name)

####################################################################################
######## Reading geometry starts here :
####################################################################################

class Read_Geo:
    def __init__(self, file_path):
#        global geometry_inserted
        self.file_path = file_path
        self.geo_frame = None
        style_colour_label = ttk.Style()
        style_colour_label.configure("Colour_Label.TLabel", foreground="black", background="Ivory", relief="flat", font=("Arial", 14))
        style_colour_frame = ttk.Style()
        style_colour_frame.configure("Colour_Frame.TFrame", background="lightblue")
        self.atoms = []  # List to store atom data
        self.symat=[]      # list to store atom names
        self.symatno=[]    # list to store atomic numbers
        self.coordx=[]     # list to store atoms x coordinates
        self.coordy=[]     # list to store atoms y coordinates
        self.coordz=[]     # list to store atoms z coordinates


        
    def Deconvert(self, value):
        # Check if the string contains 'D' for scientific notation
        if 'D' in value:
            # Replace 'D' with 'e' for Python to interpret it
            value = value.replace('D', 'e')
        
        return float(value)


    def read_geometry(self):
#        global geometry_inserted
#        Browse a file and read its lines to extract atom data.
#        Handles two formats of geometry files:
#        - 4 columns: atom name, x, y, z coordinates
#        - 5 columns: atom name, atomic number, x, y, z coordinates


        try:
            with open(self.file_path, 'r') as file:
                lines = file.readlines()


            for line in lines:
                line = line.strip()
                if not line:  # Skip empty lines
                    continue

                columns = line.split()
                if len(columns) == 4:
                    # 4-column format: atom name, x, y, z
                    atom_name = columns[0].capitalize()  # Capitalize atom name
                    x = self.Deconvert(columns[1])
                    y = self.Deconvert(columns[2])
                    z = self.Deconvert(columns[3])
#                    x, y, z = map(float, columns[1:])
                    self.atoms.append({
                        "atom": atom_name,
                        "x": x,
                        "y": y,
                        "z": z
                    })
                elif len(columns) == 5:
                    # 5-column format: atom name, atomic number, x, y, z
                    atom_name = columns[0].capitalize()  # Capitalize atom name
                    atomic_number = self.Deconvert(columns[1])
#                    x, y, z = map(float, columns[2:])
                    x = self.Deconvert(columns[2])
                    y = self.Deconvert(columns[3])
                    z = self.Deconvert(columns[4])
                    self.atoms.append({
                        "atom": atom_name,
                        "atomic_number": atomic_number,
                        "x": x,
                        "y": y,
                        "z": z
                    })
                else:
                    messagebox.showerror("Unknown Format", "Geometry file format is unknown Please Check help")
                    continue

            print('atom',self.atoms)
            total_atoms = len(self.atoms)
            for atom in self.atoms:
                self.symat.append(atom["atom"])
                self.coordx.append(atom["x"])
                self.coordy.append(atom["y"])
                self.coordz.append(atom["z"])

            for element in self.symat:
                if element in GlobVar.at_list_bold:
                    atomic_number = GlobVar.at_list_bold.index(element) + 1  # Atomic number = index + 1
                    self.symatno.append(float(atomic_number))
                else:
                    raise ValueError(f"Element {element} is not spelled Correctly.\n"
                                     " or maybe it's above 88 elements of the periodic table \n"
                                     "we only consider fistr 88 elements of the periodic table")


            print("self.symat, self.coordx, self.coordy, self.coordz, self.symatno",self.symat, self.coordx, self.coordy, self.coordz, self.symatno)

            GlobVar.geometry_inserted = True
            return self.symat, self.coordx, self.coordy, self.coordz, self.symatno

        except Exception as e:
            print(f"An error occurred: {e}")
            return None

    def insert_geo(self, root):
        geo_window = tk.Toplevel(root)
        geo_window.title("Geometry")
        geo_window.geometry("650x650")
        geo_window.configure(background="lightblue")
        entry_widgets = []


        frame1 = ttk.Frame(geo_window, style = "Colour_Frame.TFrame", padding = 10 )
        frame1.grid(row = 0, column = 0)
        frame3 = ttk.Frame(geo_window, style = "Colour_Frame.TFrame", padding = 10 )
        frame3.grid(row = 2, column = 0)
         # Scrollable container
        container = ttk.Frame(geo_window)
        container.grid(row=1, column=0, sticky=tk.W+tk.E)

        canvas = tk.Canvas(container, background="lightblue", height=450, width=600)
        scrollbar = ttk.Scrollbar(container, orient=tk.VERTICAL, command=canvas.yview)
        frame2 = ttk.Frame(canvas, style="Colour_Frame.TFrame", padding=10)

        # Configure scrollable area
        frame2.bind(
            "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=frame2, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        label = ttk.Label(frame1, text = "Please insert the number of atoms: ", style ="Colour_Label.TLabel" )
        label.grid(row = 0, column = 0, padx = 10, pady = 10, sticky = tk.W)

        atom_num_entry = ttk.Entry(frame1, width = 10)
        atom_num_entry.grid(row = 0, column = 1, padx = 10, pady = 10 )

        def create_pan():
            self.atom_num = int(atom_num_entry.get())


            for widget in frame2.winfo_children():
                widget.destroy()
            entry_widgets.clear()


            label1 = ttk.Label(frame2, text = " Atom ", background ="lightblue")
            label1.grid(row = 1, column = 0, padx = 10)
            label2 = ttk.Label(frame2, text = " X Coordinates ", background ="lightblue")
            label2.grid(row = 1, column = 1, padx = 10)
            label3 = ttk.Label(frame2, text = " Y Coordinates ", background ="lightblue")
            label3.grid(row = 1, column = 2, padx = 10)
            label4 = ttk.Label(frame2, text = " Z Coordinates ", background ="lightblue")
            label4.grid(row = 1, column = 3, padx = 10)


            for i in range (self.atom_num):
                Atom = ttk.Entry(frame2, width = 10)
                Atom.grid(row = 2+i, column = 0, padx = 10, pady = 10)

                X_coord = ttk.Entry(frame2, width = 15)
                X_coord.grid(row = 2+i, column = 1, padx = 10, pady = 10)

                Y_coord = ttk.Entry(frame2, width = 15)
                Y_coord.grid(row = 2+i, column = 2, padx = 10, pady = 10)

                Z_coord = ttk.Entry(frame2, width = 15)
                Z_coord.grid(row = 2+i, column = 3, padx = 10, pady = 10)

                entry_widgets.append((Atom, X_coord, Y_coord, Z_coord))

        def fetch_data():
#            global geometry_inserted
            for atom, x_entry, y_entry, z_entry in entry_widgets:
                atom_name = atom.get().strip()
                if atom:  # Only process if the atom field is not empty
                    x = float(self.Deconvert(x_entry.get().strip()))
                    y = float(self.Deconvert(y_entry.get().strip()))
                    z = float(self.Deconvert(z_entry.get().strip()))

                    self.atoms.append({
                        "atom": atom_name,
                        "x": x,
                        "y": y,
                        "z": z
                    })

            total_atoms = self.atom_num
            for atom in self.atoms:
               self.symat.append(atom["atom"])
               self.coordx.append(atom["x"])
               self.coordy.append(atom["y"])
               self.coordz.append(atom["z"])

            for element in self.symat:
                if element in GlobVar.at_list_bold:
                    atomic_number = GlobVar.at_list_bold.index(element) + 1  # Atomic number = index + 1
                    self.symatno.append(float(atomic_number))
                else:
                    raise ValueError(f"Element {element} is not spelled Correctly.\n"
                                     " or maybe it's above 88 elements of the periodic table \n"
                                     "we only consider fistr 88 elements of the periodic table")


            print("self.symat, self.coordx, self.coordy, self.coordz, self.symatno",self.symat, self.coordx, self.coordy, self.coordz, self.symatno)
            GlobVar.geometry_inserted = True

            return self.symat, self.coordx, self.coordy, self.coordz, self.symatno


        button = ttk.Button(frame1, text = "Enter", command = create_pan)
        button.grid(row = 0, column = 3, padx = 10, pady = 10, sticky = tk.W )

        button = ttk.Button(frame3, text = "Insert", command = fetch_data)
        button.grid(row = 0, column = 3, padx = 10, pady = 10, sticky = tk.W )
        button = ttk.Button(frame3, text = "Close", command = geo_window.destroy)
        button.grid(row = 1, column = 3, padx = 10, pady = 10, sticky = tk.W )

    def display_geometry(self):
#        global geometry_inserted
        if GlobVar.geometry_inserted == True:
            # Create a Toplevel window to display the file content
            geo_file_display = tk.Toplevel()
            geo_file_display.title("geometry File Content")
            geo_file_display.geometry("600x400")  # Optional: Set the size of the Toplevel window

            style = ttk.Style()
            style.configure("Treeview", font=("Helvetica", 14))
            # Define the Treeview widget with columns
            columns = ("atom", "x", "y", "z")
            tree = ttk.Treeview(geo_file_display, columns=columns, show="headings")
            tree.pack(expand=True, fill="both")

            # Define headings
            tree.heading("atom", text="Atom")
            tree.heading("x", text="X")
            tree.heading("y", text="Y")
            tree.heading("z", text="Z")

            # Insert rows from self.atoms
            for atom in self.atoms:
                tree.insert("", tk.END, values=(atom["atom"], atom["x"], atom["y"], atom["z"]))

            # Add a scrollbar for better usability
            scrollbar = ttk.Scrollbar(geo_file_display, orient="vertical", command=tree.yview)
            tree.configure(yscrollcommand=scrollbar.set)
            scrollbar.pack(side="right", fill="y")

    def get_geometry_data(self):
        return self.symat, self.coordx, self.coordy, self.coordz, self.symatno


###################################################################################
########## orbital section starts here :
###################################################################################

class Orb_Input:
    def __init__(self, root, keywd_button):
        self.num_orbital = GlobVar.num_orbital

        self.orbital_frame = None
        self.keywd_button = keywd_button
        self.atm_entry = []
        self.typ_entry = []
        self.atoset = np.zeros((200, 50), dtype = int)
        self.norbsym_py = np.zeros(50, dtype = int)
        self.orbsym = np.zeros((20, 20), dtype = int)
        self.activeatoms = np.zeros(30, dtype = int)
        self.atn_vector = np.zeros(200, dtype = int)
#       self.assoatm_entries = []
#       self.assotyp_entries = []
        self.orbital_data = []
        self.description_orb = ""
        self.create_orbital_section()
                
#        style_colour_label1 = ttk.Style()
#        style_colour_label1.configure("Colour_Label1.TLabel", foreground="black", background="lightblue", relief="flat", font=("Arial", 16))


                # Insert button to validate and save inputs
        insert_button = ttk.Button(
            self.orbital_frame, text="Insert", command=self.validate_and_store_orbital_data
        )
        insert_button.grid(row=3, column=0, columnspan=2, pady=10)
        close_button = ttk.Button(self.orbital_frame, text = "Close", command = self.orbital_frame.destroy)
        close_button.grid(row = 4, column = 0, pady = 10, columnspan=2 )

    def create_orbital_section(self):
        # Create a new frame for orbital inputs if it doesn't exist
        if self.orbital_frame is None:
            self.orbital_frame = tk.Toplevel(root, padx=10, pady=10)
            self.orbital_frame.title("orbital inputs")
            self.orbital_frame.geometry("560x560")
            self.orbital_frame.configure(background="lightblue")

        for widget in self.orbital_frame.winfo_children():
            widget.destroy()

        frame = ttk.Frame(self.orbital_frame, style = "Colour_Frame.TFrame")
        frame.grid(row = 0, column = 0)
        frame0 = ttk.Frame(self.orbital_frame, style = "Colour_Frame.TFrame")
        frame0.grid(row = 1, column = 0)

         # Scrollable container
        container = ttk.Frame(self.orbital_frame)
        container.grid(row=2, column=0, sticky=tk.W+tk.E)

        canvas = tk.Canvas(container, background="lightblue", height=370, width=510)
        scrollbar = ttk.Scrollbar(container, orient=tk.VERTICAL, command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)

        frame1 = ttk.Frame(canvas, style = "Colour_Frame.TFrame")
        canvas.create_window((0, 0), window=frame1, anchor="nw")

        # Configure scrollable area
        frame1.bind(
            "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=frame1, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Create input panes for each active orbitals
        ttk.Label(frame, text=f"Number of active orbitals: {self.num_orbital}", style = "Colour_Label.TLabel" ).grid(row=0, column=0, columnspan=3, pady=15)

        #creating the pane for getting orbital numbers
        ttk.Label(frame0, text=f"     Atom Number", style = "Colour_Label1.TLabel").grid(row=0, column=3, padx=3, sticky=tk.E)
        ttk.Label(frame0, text=f"Orbital type", style = "Colour_Label1.TLabel").grid(row=0, column=4, padx=5,  sticky=tk.E)

#        ttk.Label(frame1, text=f"--------------------", style = "Colour_Label1.TLabel").grid(row=2, column=1, padx=5, pady=5, sticky=tk.W)
#        ttk.Label(frame1, text=f"--------------------", style = "Colour_Label1.TLabel").grid(row=2, column=2, padx=5, pady=5, sticky=tk.W)

        for i in range(self.num_orbital):
            ttk.Label(frame1, text=f"active orbital {i+1} ", style = "Colour_Label.TLabel" ).grid(row=i+3, column=0, padx=30, pady=10, sticky=tk.W)
#       creating the pane for getting associated atom number
            self.assoatm_entry = ttk.Entry(frame1, width=10)
            self.assoatm_entry.grid(row=i+3, column=1, padx=30, pady=10)
            self.atm_entry.append(self.assoatm_entry)

#       creating the pane for getting associated atom number
            self.assotyp_entry = ttk.Entry(frame1, width=10)
            self.assotyp_entry.grid(row=i+3, column=2, padx=30, pady=10, sticky = tk.E)
            self.typ_entry.append(self.assotyp_entry)


    def validate_and_store_orbital_data(self):
        """Validate orbital inputs and store them in a list."""
        self.orbital_data.clear()  # Clear previous data

        for i in range(self.num_orbital):
            atom_number = self.atm_entry[i].get().strip()
            orbital_type = self.typ_entry[i].get().strip()

            # Validate inputs
            if not atom_number or not orbital_type:
                messagebox.showerror(
                    "Validation Error",
                    f"Please fill in all fields for orbital {i + 1}.",
                )
                self.atm_entry[i].configure(background="pink")
                self.typ_entry[i].configure(background="pink")
                return

            try:
                atom_number = int(atom_number)  # Ensure atom number is an integer
            except ValueError:
                messagebox.showerror(
                    "Input Error", f"Atom number for orbital {i + 1} must be an integer."
                )
                self.atm_entry[i].configure(background="pink")
                return

            # Reset background color after validation
            self.atm_entry[i].configure(background="white")
            self.typ_entry[i].configure(background="white")

            # Store the data
            self.orbital_data.append({
                "atom_number": atom_number,
                "orbital_type": orbital_type
            })

        messagebox.showinfo("Success", "Orbital inputs validated and stored successfully.")
        sorbs = []
        pxorbs = []
        pyorbs = []
        pzorbs = []
        i = 0
        sig_type = 0
        px_type = 0
        py_type = 0
        pz_type = 0
        for data in self.orbital_data:
            i=i+1
            if 's' in data["orbital_type"].lower():
                sig_type = sig_type+1
                sorbs.append(i+GlobVar.num_iao)
            elif 'px' in data["orbital_type"].lower():
                px_type = px_type + 1
                pxorbs.append(i+GlobVar.num_iao)
            elif 'py' in data["orbital_type"].lower():
                py_type = py_type + 1
                pyorbs.append(i+GlobVar.num_iao)
            elif 'pz' in data["orbital_type"].lower():
                pz_type = pz_type + 1
                pzorbs.append(i+GlobVar.num_iao)
            elif 'pi1' in data["orbital_type"].lower():
                px_type = px_type + 1
                pxorbs.append(i+GlobVar.num_iao)
            elif 'pi2' in data["orbital_type"].lower():
                py_type = py_type + 1
                pyorbs.append(i+GlobVar.num_iao)
            elif 'pi3' in data["orbital_type"].lower():
                pz_type = pz_type + 1
                pzorbs.append(i+GlobVar.num_iao)
            else:
                messagebox.showerror("Input Error",f"The type of the orbital {i} is unknown. \n"
                                     "Please put s, px, py or pz otherwise put sig or sigma,\n"
                                     " pi1, pi2, pi3 if the direction of pi orbs are different"
                        )
        norbsym = [
                px_type, 
                py_type, 
                pz_type, 
                sig_type
                ]

        self.norbsym_py = np.zeros(50, dtype=int)
        self.norbsym_py[:len(norbsym)]=norbsym
#        max_len = max(len(pxorbs), len(pyorbs), len(pzorbs), len(sorbs))
        max_len = 20 

        # Pad lists to be of equal length
        pxorbs_padded = pxorbs + [0] * (max_len - len(pxorbs))
        pyorbs_padded = pyorbs + [0] * (max_len - len(pyorbs))
        pzorbs_padded = pzorbs + [0] * (max_len - len(pzorbs))
        sorbs_padded = sorbs + [0] * (max_len - len(sorbs))
        self.orbsym[0,:] = pxorbs_padded
        self.orbsym[1,:] = pyorbs_padded
        self.orbsym[2,:] = pzorbs_padded
        self.orbsym[3,:] = sorbs_padded

        print('norbsym', self.norbsym_py)
        print('sig_type, px_type, py_type, pz_type',sig_type, px_type, py_type, pz_type )
        # Check how many types are non-zero
        GlobVar.type_orb_count = sum(1 for count in [sig_type, px_type, py_type, pz_type] if count > 0)

        self.atoset = self.create_matrix() 
        "atoset is matrix where each row represents an atom according to the geometry, if the atom ia an active atom,"
        "the corresponding column get '1' otherwise '0'. and the next columns contain the corresponding active orbital"
        "numbers associated with that atom."
        print('atoset',self.atoset)
        self.keywd_button.config(state=tk.NORMAL)  

    def create_matrix(self):
        # Step 1: Gather atom numbers and their orbitals
        print('self.orbital_data',self.orbital_data)
        atom_to_orbitals = {}
        active=[]
        for orbital_number, entry in enumerate(self.orbital_data, start=1):
            atom_number = entry["atom_number"] 
            if atom_number not in atom_to_orbitals:
                atom_to_orbitals[atom_number] = []
            atom_to_orbitals[atom_number].append(orbital_number)
        print('atom_to_orbitals',atom_to_orbitals)

        # Step 2: Determine the maximum atom number
        max_atom_number = max(atom_to_orbitals.keys())
        max_orbitals = max(len(orbitals) for orbitals in atom_to_orbitals.values())

        matrix = np.zeros((200, 20), dtype=int)

        for atom_number, orbitals in atom_to_orbitals.items():
            matrix[atom_number - 1, 0] = 1  # Mark presence
            for col_index, orbital in enumerate(orbitals):
                matrix[atom_number - 1, col_index + 1] = orbital + GlobVar.num_iao # num_iao = number of inactive orbitals

        active = np.where(matrix[:, 0] == 1)[0]+1

        self.activeatoms[:len(active)]= active
        print(active, self.activeatoms)

        atn_vec = [sum(1 for x in row if x != 0) for row in matrix]
        self.atn_vector[:len(atn_vec)]=atn_vec

        return matrix

    def get_orbital_matrices(self):
        atoset_matrix = self.atoset
        norbsym_vector= self.norbsym_py
        active_atoms = self.activeatoms
        atn = self.atn_vector
        orbsym_matrix = self.orbsym
        return (atoset_matrix, norbsym_vector, active_atoms, atn, orbsym_matrix)


class Keywd_Input:
    def __init__(self, root, Run_Button):
        self.root = root
        self.Run_Button = Run_Button
        self.multiplicity = GlobVar.multiplicity
        self.type_orb_count = GlobVar.type_orb_count
        self.structure_type_entry = False
        self.rumer_set_num_entry = False
        self.method_type_entry = False
        self.keywd_window = None
        self.priority_window = None
        self.priority_str_window = None
        self.prio_bond_window = None
        self.prio_rads_window = None
        self.frame = None
        self.method_type = tk.StringVar(value="Chem inst")
        self.ChemInst_set_type = tk.StringVar(value="Single Set")
        self.rumer_set_type = tk.StringVar(value="Single Rumer Set")
        self.str_type = tk.StringVar(value="Covalent")
        self.cheminst_str_type = tk.StringVar(value="Symmetry")
        self.symmetry_set_order_type = tk.StringVar(value="Quality-Arrange")
        self.ovlp_cal = tk.StringVar(value="No")
        self.IAB_type = tk.StringVar(value="1")
        self.NAB_type = tk.StringVar(value="2")
        self.SBB_type = tk.StringVar(value="3")
        self.PDB_type = tk.StringVar(value="None")
        self.PDR_type = tk.StringVar(value="None")
        self.ChemInst_set_type_entry = True
        self.cheminst_str_type_entry = False
        self.symmetric_set_order_type_entry = False
        self.prio_structure_entries = []  
        self.prio_bond_entries = []  
        self.prio_rads_entries = []  
        self.PDR_buttons = []
        self.SBB_buttons = []
        self.PDB_data = []
        Pmw.initialise(self.root)
        self.balloon = Pmw.Balloon(self.root)
        self.mout_number = 1
        self.bond_number = 0
        self.Run_Button.config(state=tk.NORMAL)  

    def create_keywd_pane(self):
        if self.keywd_window is None:
            self.keywd_window = tk.Toplevel(self.root, padx = 10, pady = 10)
            self.keywd_window.title("Spatial Keyword Inputs")
            self.keywd_window.geometry("750x850")
            self.keywd_window.configure( background = "lightblue")

            self.method_type_frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.method_type_frame.grid (row = 0, column = 0, sticky = tk.W)

            self.rumer_set_type_frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.rumer_set_type_frame.grid (row = 1, column = 0, sticky = tk.W)

            self.str_type_frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.str_type_frame.grid (row = 2, column = 0, sticky = tk.W)

#            self.non_rum_frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
#            self.non_rum_frame.grid (row = 3, column = 0, sticky = tk.W)

            self.ChemInst_set_type_frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.ChemInst_set_type_frame.grid(row = 3, column = 0, sticky = tk.W)

            self.cheminst_str_type_frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.cheminst_str_type_frame.grid (row = 4, column = 0, sticky = tk.W)

            self.symmetry_set_order_frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.symmetry_set_order_frame.grid (row = 5, column = 0, sticky = tk.W)

            self.ovlp_cal_frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.ovlp_cal_frame.grid (row = 6, column = 0, sticky = tk.W)

            self.prio_label_frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.prio_label_frame.grid (row = 7, column = 0, sticky = tk.W)

            self.diff_qualities_frame = ttk.Frame(self.keywd_window, width = 30,  style="Colour_Frame.TFrame", padding = 10)
            self.diff_qualities_frame.grid(row = 8, column=0, sticky="nsew", columnspan = 2)

            self.PDB_PDR_Button_frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.PDB_PDR_Button_frame.grid (row = 9, column = 0, sticky = tk.W)

            self.priority_str_frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.priority_str_frame.grid(row = 10, column = 0, sticky = tk.W)

            self.mout_frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.mout_frame.grid (row = 11, column = 0, sticky = tk.W)

            self.end_frame = ttk.Frame(self.keywd_window, width = 30,  style="Colour_Frame.TFrame", padding = 10)
            self.end_frame.grid(row=12, column=0, sticky="nsew", columnspan = 2)

#            self.frame9 = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
#            self.frame9.grid (row = 10, column = 0)

            method_type_label = ttk.Label(self.method_type_frame, text="Select method type: ", style = "Colour_Label.TLabel")
            method_type_label.grid(row = 0, column = 0, sticky = tk.W, padx = 10, pady=10)
            self.balloon.bind(method_type_label,"Select any one from the two method provided. \n"
                                                    "select 'Rumer' to genarate a Rumer structure set,\n"
                                                    " no need to provide any spatial keywords. For all\n" 
                                                    "other type of set please select Chem Inst.")

            Rumer_set_type_label = ttk.Label(self.rumer_set_type_frame, text="Rumer Set type: ", style = "Colour_Label.TLabel")
            Rumer_set_type_label.grid(row = 0, column = 0, sticky = tk.W, padx = 10, pady=10)
            self.balloon.bind(Rumer_set_type_label,"Select any one from the two Rumer set options. \n"
                                        "select 'Single Rumer Set' to genarate one Rumer structure set,\n"
                                        " corresponding to the given order of orbitals. Select 'All Rumer Sets'\n" 
                                        "to get all possible unique Rumer sets.")

            str_type_label = ttk.Label(self.str_type_frame, text="Select structure type: ", style = "Colour_Label.TLabel")
            str_type_label.grid(row = 0, column = 0, sticky = tk.W, padx = 10, pady=10)
            self.balloon.bind(str_type_label,"Select any one from the three types. select 'Covalent'\n"
                                         " to generate covalent sets only. Select 'Ionic' to generate \n" 
                                         "ionic structures sets only. Select both to generate both \n"
                                          "type of structure sets.")

            ChemInst_set_type_label = ttk.Label(self.ChemInst_set_type_frame, text = "Chem Inst Set type:", style = "Colour_Label.TLabel")
            ChemInst_set_type_label.grid(row = 0, column = 0, sticky = tk.W, padx = 10, pady = 10)
            self.balloon.bind(ChemInst_set_type_label,"Select any one from the four types of Chemical insight output set.\n"
                              " select 'Single Set' if you want the best chemically insightfull one set,\n" 
                              " select 'All Best Sets' if you want all possible same best quality sets if available,\n" 
                              " select 'All Sets' if you want all possible sets but it can be huge adjust maximum \n" 
                              " output file accordingly, select 'Eq Bond' if you want only equally distributed bond sets\n"
                              " in other words this sets are lowest overlapped sets but less chemically meaningfull. \n" 
                              " these sets could be helpfull to read off negetive Coulson-Chirgwin weights")

            ChemInst_str_type_label = ttk.Label(self.cheminst_str_type_frame, text = "Cheminst Str Type:", style = "Colour_Label.TLabel")
            ChemInst_str_type_label.grid(row = 0, column = 0, padx = 10, pady = 10)
            self.balloon.bind(ChemInst_str_type_label,"Select any one from the three types of calculations. \n"
                              "Select 'Symmetry' if you want to have symmetric sets. \n"
                              "Select 'Asymetric' if you want to have chemical insight sets. \n"
                              "Select 'Checksymm' if you want to check some sets are symmetric or not")

            symmetry_set_order_label = ttk.Label(self.symmetry_set_order_frame, text = "Symmetric group arrangement:", style = "Colour_Label.TLabel")
            symmetry_set_order_label.grid(row = 0, column = 0, sticky = tk.W, padx = 10, pady = 10)
            self.balloon.bind(symmetry_set_order_label,"Select any order of arrangement of the symmetric group. \n"
                              "It could help searching a symmetric set easyer. 'Quality' arrange"
                              "symmetric groups from higher quality to lower quality. other two options\n"
                              "'big to small' and 'small to big' arrange the groups according to their sizes")

            ovlp_cal_label = ttk.Label(self.ovlp_cal_frame, text = "Estimate Overlap", style = "Colour_Label.TLabel")
            ovlp_cal_label.grid(row=0,column = 0, padx=10, pady=10, columnspan = 2)
            self.balloon.bind(ovlp_cal_label, "To estimate overlap among the structures \n"
                                          "in each set click yes. The default is set to 'No'")

            tip_method = ["Chem inst", "Rumer"]

            for i, tip in enumerate(tip_method, start=1):
                button = ttk.Radiobutton(
                    self.method_type_frame,
                    text=tip,
                    value=tip,
                    variable=self.method_type,
                    command=self.update_method_type,
                    style="Custom.TRadiobutton"
                )
                button.grid(row = 0, column = i, padx=10, pady=10)

            Rumer_Set_Type = ["Single Rumer Set", "All Rumer Sets"]

            for i, Set in enumerate(Rumer_Set_Type, start=1):
                button = ttk.Radiobutton(
                    self.rumer_set_type_frame,
                    text=Set,
                    value=Set,
                    variable=self.rumer_set_type,
                    command=self.Update_Rumer_Set_Type,
                    style="Custom.TRadiobutton"
                )
                button.grid(row = 0, column = i, padx=10, pady=10)
                self.disable_widgets_in_frame(self.rumer_set_type_frame)

            str_type = ["Covalent", "Ionic", "Both"]
            # Create and pack each radio button
            for i, tip in enumerate(str_type, start=1):
                button = ttk.Radiobutton(
                    self.str_type_frame,
                    text=tip,
                    value=tip,
                    variable=self.str_type,
                    command=self.Update_Str_Type,
                    style="Custom.TRadiobutton"
                )
                button.grid(row = 0, column = i, padx=10, pady=10)

            cheminstsets_type = ["Single Set", "All Best Sets", "All Sets", "Eq Bond"]
            for i, sets in enumerate(cheminstsets_type, start=1):
                button2 = ttk.Radiobutton(
                    self.ChemInst_set_type_frame,
                    text=sets,
                    value=sets,
                    variable=self.ChemInst_set_type,
                    command=self.ChemInst_set_type_read,
                    style="Custom.TRadiobutton",
                )
                button2.grid(row=0, column=i, padx=10, pady=10)


            Cheminst_str_type = ["Symmetry", "Asymmetry", "Checksymm"]
            for i, chem in enumerate(Cheminst_str_type, start=1):
                button1 = ttk.Radiobutton(
                    self.cheminst_str_type_frame,
                    text=chem,
                    value=chem,
                    variable=self.cheminst_str_type,
                    command=self.cheminst_str_type_read,
                    style="Custom.TRadiobutton",
                )
                button1.grid(row=0, column=i, padx=10, pady=10)

            set_order = ["Quality-Arrange","Big-to-Small","Small-to-Big"]
            for i, seto in enumerate(set_order, start=1):
                button3 = ttk.Radiobutton(
                    self.symmetry_set_order_frame,
                    text=seto,
                    value=seto,
                    variable=self.symmetry_set_order_type,
                    command=self.symmetry_set_order_type_read,
                    style="Custom.TRadiobutton",
                )
                button3.grid(row=0, column=i, padx=10, pady=10)



            ovlp_types = ["Yes", "No"]
            for i, ovlp in enumerate(ovlp_types, start=2):
                ovlp_button = ttk.Radiobutton(
                    self.ovlp_cal_frame,
                    text=ovlp,
                    value=ovlp,
                    variable=self.ovlp_cal,
                    command=self.get_str_ovlp_keywd,
                    style="Custom.TRadiobutton",
                )
                ovlp_button.grid(row=0, column=i, padx=10, pady=10)

            prio_label = ttk.Label(self.prio_label_frame, text = " Decide Priorities of the Chemical Qualities: ", style = "Colour_Label.TLabel")
            prio_label.grid(row = 0, column = 0, padx = 10, pady = 5)
            self.balloon.bind(prio_label, "The sequence of the default priority is IAB > NAB > SBB and \n"
                              "PBU & PRU are not taken in the calculation; To change this default priorities\n"
                              "please click the button and selevct the numbers for each ")

            IAB_label = ttk.Label(self.diff_qualities_frame, text = " IAB  PRIORITY: ", style = "Colour_Label2.TLabel")
            IAB_label.grid(row = 0, column = 0, padx = 10, pady = 5)
            self.balloon.bind(IAB_label, "Intra Atomic Bond Priority. Default is set to '1'")

            IAB_option = ["1", "2", "3", "4", "5", "None"]
            for i, IAB in enumerate(IAB_option, start=1):
                IAB_button = ttk.Radiobutton(
                    self.diff_qualities_frame,
                    text=IAB,
                    value=IAB,
                    variable=self.IAB_type,
                    command=self.get_IAB_Priority,
                    style="Custom.TRadiobutton",
                )
                IAB_button.grid(row=0, column=i, padx=10, pady=10)

            NAB_label = ttk.Button(self.diff_qualities_frame, text = " NAB PRIORITY: ", style = "Colour_Label2.TLabel")
            NAB_label.grid(row = 1, column = 0, padx = 10, pady = 5)
            self.balloon.bind(NAB_label, "Near Atomic Bond Priority. Default is set to '2'")

            NAB_option = ["1", "2", "3", "4", "5", "None"]
            for i, NAB in enumerate(NAB_option, start=1):
                NAB_button = ttk.Radiobutton(
                    self.diff_qualities_frame,
                    text=NAB,
                    value=NAB,
                    variable=self.NAB_type,
                    command=self.get_NAB_Priority,
                    style="Custom.TRadiobutton",
                )
                NAB_button.grid(row=1, column=i, padx=10, pady=10)

            SBB_label = ttk.Button(self.diff_qualities_frame, text = " SBB PRIORITY: ", style = "Colour_Label2.TLabel")
            SBB_label.grid(row = 2, column = 0, padx = 10, pady = 5)
            self.balloon.bind(SBB_label, "Symmetry Breaking Bond Priority. Default is set to '3'")

            SBB_option = ["1", "2", "3", "4", "5", "None"]
            for i, SBB in enumerate(SBB_option, start=1):
                SBB_button = ttk.Radiobutton(
                    self.diff_qualities_frame,
                    text=SBB,
                    value=SBB,
                    variable=self.SBB_type,
                    command=self.get_SBB_Priority,
                    style="Custom.TRadiobutton",
                )
                SBB_button.grid(row=2, column=i, padx=10, pady=10)
                self.SBB_buttons.append(SBB_button)

            PDB_label = ttk.Button(self.diff_qualities_frame, text = " PDB PRIORITY: ", style = "Colour_Label2.TLabel")
            PDB_label.grid(row = 3, column = 0, padx = 10, pady = 5)
            self.balloon.bind(PDB_label, "Pre-defined Bond Priority. Default is set to \n"
                              "'None': Not taken into the calculation")

            PDB_option = ["1", "2", "3", "4", "5", "None"]
            for i, PDB in enumerate(PDB_option, start=1):
                PDB_button = ttk.Radiobutton(
                    self.diff_qualities_frame,
                    text=PDB,
                    value=PDB,
                    variable=self.PDB_type,
                    command=self.get_PDB_Priority,
                    style="Custom.TRadiobutton",
                )
                PDB_button.grid(row=3, column=i, padx=10, pady=10)

            PDR_label = ttk.Button(self.diff_qualities_frame, text = " PDR PRIORITY: ", style = "Colour_Label2.TLabel")
            PDR_label.grid(row = 4, column = 0, padx = 10, pady = 5)
            self.balloon.bind(PDR_label, "Pre-defined Bond Priority. Default is set to \n"
                              "'None': Not taken into the calculation")

            PDR_option = ["1", "2", "3", "4", "5", "None"]
            for i, PDR in enumerate(PDR_option, start=1):
                PDR_button = ttk.Radiobutton(
                    self.diff_qualities_frame,
                    text=PDR,
                    value=PDR,
                    variable=self.PDR_type,
                    command=self.get_PDR_Priority,
                    style="Custom.TRadiobutton",
                )
                PDR_button.grid(row=4, column=i, padx=10, pady=10)

                self.PDR_buttons.append(PDR_button)
            self.check_and_disable_buttons()


            self.PDB_button = ttk.Button(self.PDB_PDR_Button_frame, text = " Insert PDB ", command = self.Insert_PDB)
            self.PDB_button.grid(row = 0, column = 0, padx = 10, sticky = tk.W )
            self.PDB_button.config(state=tk.DISABLED)  # Initially disable the button

            self.PDR_button = ttk.Button(self.PDB_PDR_Button_frame, text = " Insert PDR ", command = self.Insert_PDR)
            self.PDR_button.grid(row = 0, column = 1, padx = 10, sticky = tk.W )
            self.PDR_button.config(state=tk.DISABLED)  # Initially disable the button

            prio_struc_button = ttk.Button(self.priority_str_frame, text = "Insert Priorities Structures", command = self.Insert_Priority_Str )
            prio_struc_button.grid(row = 0, column = 0, padx = 10, pady = 5)
            self.balloon.bind(prio_struc_button, "If you want some structures must be present in the set \n"
                                                 "you can insert them by clicking this Button;\n" 
                                                 "The provided structures must be linearly independent \n"
                                                 "otherwise they will not be present in the set together")
            
            mout_label = ttk.Label(self.mout_frame, text = "Maximum Number of Output Files:", style = "Colour_Label.TLabel")
            mout_label.grid(row=1,column = 0, padx=10, pady=10, columnspan = 2)
            mout_entry = ttk.Entry(self.mout_frame, width = 10)
            mout_entry.grid(row = 1, column= 2, padx = 10, pady = 10)
            mout_entry.insert(0, "1")
            self.balloon.bind(mout_label,"For a large Active Space, the total number of sets can\n"
                                         " reach millions or more. Each output file can contain\n" 
                                         " up to 75,000 sets. By default, the number of output files\n"
                                         " is set to one. If additional output files are required, please\n"
                                         " specify the desired number in the entry box and press 'Enter'.")
            

            insert_button = ttk.Button(self.end_frame, text = "Insert", command = self.get_keywds)
            insert_button.grid(row = 0, column = 0, pady = 10, columnspan=2, sticky = tk.W )

            close_button = ttk.Button(self.end_frame, text = "Close", command = self.keywd_window.destroy)
            close_button.grid(row = 1, column = 0, pady = 10, columnspan=2, sticky = tk.W )

    def check_and_disable_buttons(self):
 # Disable PDR options if multiplicity = 1, means the system is singlet
        print('multiplicity',self.multiplicity)
        if self.multiplicity == 1:
            for button in self.PDR_buttons:
                button.config(state=tk.DISABLED)   
        if self.type_orb_count == 1:
            for button in self.SBB_buttons:
                button.config(state=tk.DISABLED)   

    def get_IAB_Priority(self):
        IAB_Priority = self.IAB_type.get()
        print('IAB_Priority', IAB_Priority)
        return (IAB_Priority)

    def get_NAB_Priority(self):
        NAB_Priority = self.NAB_type.get()
        print('NAB_Priority', NAB_Priority)
        return (NAB_Priority)

    def get_SBB_Priority(self):
        SBB_Priority = self.SBB_type.get()
        print('SBB_Priority', SBB_Priority)
        return (SBB_Priority)

    def get_PDB_Priority(self):
        PDB_Priority = self.PDB_type.get()
        print('PDB_Priority', PDB_Priority)
        if PDB_Priority is not None:
            self.PDB_button.config(state=tk.NORMAL)
        return (PDB_Priority)

    def get_PDR_Priority(self):
        PDR_Priority = self.PDR_type.get()
        print('PDR_Priority', PDR_Priority)
        if PDR_Priority is not None and self.multiplicity != 1:
            self.PDR_button.config(state=tk.NORMAL)
        return (PDR_Priority)

    def update_method_type(self):
        method_type = self.method_type.get()
        if method_type:
            self.method_type_entry = True
            if method_type == 'Rumer':
                self.enable_widgets_in_frame(self.rumer_set_type_frame)
                self.disable_widgets_in_frame(self.ChemInst_set_type_frame) 
                self.disable_widgets_in_frame(self.cheminst_str_type_frame)
                self.disable_widgets_in_frame(self.symmetry_set_order_frame)
                self.disable_widgets_in_frame(self.ovlp_cal_frame)
                self.disable_widgets_in_frame(self.prio_label_frame) 
                self.disable_widgets_in_frame(self.diff_qualities_frame)
                self.disable_widgets_in_frame(self.PDB_PDR_Button_frame)
                self.disable_widgets_in_frame(self.priority_str_frame)
                self.disable_widgets_in_frame(self.mout_frame)
            elif method_type == 'Chem inst':
                self.disable_widgets_in_frame(self.rumer_set_type_frame)
                self.enable_widgets_in_frame(self.ChemInst_set_type_frame) 
                self.enable_widgets_in_frame(self.cheminst_str_type_frame)
                self.enable_widgets_in_frame(self.symmetry_set_order_frame)
                self.enable_widgets_in_frame(self.ovlp_cal_frame)
                self.enable_widgets_in_frame(self.prio_label_frame) 
                self.enable_widgets_in_frame(self.diff_qualities_frame)
                self.enable_widgets_in_frame(self.PDB_PDR_Button_frame)
                self.enable_widgets_in_frame(self.priority_str_frame)
                self.enable_widgets_in_frame(self.mout_frame)
            print('method_type',method_type)
            return (method_type)

    def disable_widgets_in_frame(self,frame):
        """
        Recursively disable all widgets and frames inside the given frame.
        """
        for widget in frame.winfo_children():
            widget.configure(state="disabled")

    def enable_widgets_in_frame(self,frame):
        """
        Recursively enable all widgets and frames inside the given frame.
        """
        for widget in frame.winfo_children():
            widget.configure(state="normal")

    def Update_Rumer_Set_Type(self):
        rumer_set_num = self.rumer_set_type.get()
        if rumer_set_num:
            self.rumer_set_num_entry = True
            print('rumer_set_num',rumer_set_num)
            return (rumer_set_num)

    def Update_Str_Type(self):
        structure_type = self.str_type.get()
        if structure_type:
            self.structure_type_entry = True
            print('structure_ype',structure_type)
            return (structure_type)


    def ChemInst_set_type_read(self):
        set_type = self.ChemInst_set_type.get()
        if set_type:
            self.ChemInst_set_type_entry = True
            print('set_type',set_type)
           # if set_type == 'All Sets' or set_type =='Eq Bond' or set_type =='All Best Sets':

           #     ovlp_close_button = ttk.Button(self.frame5, text = "Close", command = self.frame5.destroy)
           #     ovlp_close_button.grid(row = 0, column= 4, padx = 10, pady = 10)
            return (set_type)

    def get_maximum_num_output(self, mout_entry):
        try:
            self.mout_number= int(mout_entry.get())
            return(self.mout_number)
        except ValueError:
            # Handle invalid input gracefully
            tk.messagebox.showerror("Invalid Input", "Please enter a valid number.")
            return
        print('mout_number', self.mout_number)
#        self.frame4.destroy()

    def get_str_ovlp_keywd(self):
        ovlp_calc = self.ovlp_cal.get()
        print('ovlp_cal', ovlp_calc)
        return (ovlp_calc)


    def cheminst_str_type_read(self):
        cheminst_type = self.cheminst_str_type.get()
        if cheminst_type:
            self.cheminst_str_type_entry = True
            print('cheminst_type',cheminst_type)
        return (cheminst_type)
        #    return (cheminst_type)
#            if cheminst_type == 'Symmetry':


    def symmetry_set_order_type_read(self):
        set_order_type = self.symmetry_set_order_type.get()
        if set_order_type:
            self.symmetry_set_order_type_entry = True
            print('set_order_type',set_order_type)
            return (set_order_type)

#    def create_priority_keywd(self):

    def Insert_PDB(self):
        if self.prio_bond_window is None:
            self.prio_bond_window = tk.Toplevel(self.root, padx=10, pady=10)
            self.prio_bond_window.title("Priority Bonds")
            self.prio_bond_window.geometry("450x450")
            self.prio_bond_window.configure(background="lightblue")

            frame = ttk.Frame(self.prio_bond_window, style="Colour_Frame.TFrame")
            frame.grid(row=0, column=0, sticky=tk.W)

             # Scrollable frame for dynamic fields
            canvas = tk.Canvas(self.prio_bond_window, background="lightblue", height=300, width=300)
            scrollbar = ttk.Scrollbar(self.prio_bond_window, orient=tk.VERTICAL, command=canvas.yview)
            scrollable_frame = ttk.Frame(canvas, style = "Colour_Frame.TFrame")

            # Configure scroll region
            scrollable_frame.bind(
                "<Configure>",
                lambda e: canvas.configure(scrollregion=canvas.bbox("all")),
            )
            canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
            canvas.configure(yscrollcommand=scrollbar.set)

            # Place the scrollable frame and scrollbar
            canvas.grid(row=1, column=0, sticky=tk.W + tk.E)
            scrollbar.grid(row=1, column=1, sticky=tk.N + tk.S)

            frame1 = ttk.Frame(self.prio_bond_window, style="Colour_Frame.TFrame")
            frame1.grid(row=2, column=0)

            label = ttk.Label(frame, text = "Number of Structures:", style = "Colour_Label.TLabel")
            label.grid(row = 0, column = 0 , padx = 10, pady = 10)
            prio_bond_entry = ttk.Entry(frame, width = 10 )
            prio_bond_entry.grid(row = 0, column = 1)

            Insert_button = ttk.Button(frame, text = "Insert", command =lambda:self.generate_prio_bond_fields(scrollable_frame, prio_bond_entry))
            Insert_button.grid(row = 0, column = 2, padx = 10, pady = 10 )

            close_button = ttk.Button(frame1, text = "DONE", command = self.get_prio_bond_data)
            close_button.grid(row = 0, column = 0, pady = 10, columnspan=2 )

    def generate_prio_bond_fields(self, frame1, prio_bond_entry):
        """
        Generate labels and entry fields dynamically based on user input.
        """
        try:
            self.bond_number = int(prio_bond_entry.get())  # Get the number of structures
        except ValueError:
            # Handle invalid input gracefully
            tk.messagebox.showerror("Invalid Input", "Please enter a valid number.")
            return

        # Clear previous fields
        for widget in frame1.winfo_children():
            widget.destroy()

        for i in range(self.bond_number):
            prio_bond_label = ttk.Label(frame1, text=f"Pre-defined Bond {i+1}:", style="Colour_Label1.TLabel")
            prio_bond_label.grid(row=i, column=0, padx=10, pady=10)

            prio_bond_entry = ttk.Entry(frame1, width=20)
            prio_bond_entry.grid(row=i, column=1, padx=10, pady=10)

            self.prio_bond_entries.append(prio_bond_entry)  # Save reference

    def get_prio_bond_data(self):
        """
        Retrieve data from the dynamically created structure entry fields.
        """
        self.PDB_data = [entry.get() for entry in self.prio_bond_entries]
        print("Entered prio_bond:", self.PDB_data)
        self.prio_bond_window.destroy()
        self.prio_bond_window = None
        return self.PDB_data

    def Insert_PDR(self):
        if self.prio_rads_window is None:
            self.prio_rads_window = tk.Toplevel(self.root, padx=10, pady=10)
            self.prio_rads_window.title("Priority Radicals")
            self.prio_rads_window.geometry("450x450")
            self.prio_rads_window.configure(background="lightblue")

            frame = ttk.Frame(self.prio_rads_window, style="Colour_Frame.TFrame")
            frame.grid(row=0, column=0, sticky=tk.W)

             # Scrollable frame for dynamic fields
            canvas = tk.Canvas(self.prio_rads_window, background="lightblue", height=300, width=300)
            scrollbar = ttk.Scrollbar(self.prio_rads_window, orient=tk.VERTICAL, command=canvas.yview)
            scrollable_frame = ttk.Frame(canvas, style = "Colour_Frame.TFrame")

            # Configure scroll region
            scrollable_frame.bind(
                "<Configure>",
                lambda e: canvas.configure(scrollregion=canvas.bbox("all")),
            )
            canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
            canvas.configure(yscrollcommand=scrollbar.set)

            # Place the scrollable frame and scrollbar
            canvas.grid(row=1, column=0, sticky=tk.W + tk.E)
            scrollbar.grid(row=1, column=1, sticky=tk.N + tk.S)

            frame1 = ttk.Frame(self.prio_rads_window, style="Colour_Frame.TFrame")
            frame1.grid(row=2, column=0)

            label = ttk.Label(frame, text = "Number of radicals:", style = "Colour_Label.TLabel")
            label.grid(row = 0, column = 0 , padx = 10, pady = 10)
            prio_rads_entry = ttk.Entry(frame, width = 10 )
            prio_rads_entry.grid(row = 0, column = 1)

            Insert_button = ttk.Button(frame, text = "Insert", command =lambda:self.generate_prio_rads_fields(scrollable_frame, prio_rads_entry))
            Insert_button.grid(row = 0, column = 2, padx = 10, pady = 10 )

            close_button = ttk.Button(frame1, text = "DONE", command = self.get_prio_rads_data)
            close_button.grid(row = 0, column = 0, pady = 10, columnspan=2 )

    def generate_prio_rads_fields(self, frame1, prio_rads_entry):
        """
        Generate labels and entry fields dynamically based on user input.
        """
        try:
            rads_number = int(prio_rads_entry.get())  # Get the number of structures
        except ValueError:
            # Handle invalid input gracefully
            tk.messagebox.showerror("Invalid Input", "Please enter a valid number.")
            return

        # Clear previous fields
        for widget in frame1.winfo_children():
            widget.destroy()

        for i in range(rads_number):
            prio_rads_label = ttk.Label(frame1, text=f"Pre-defined rad {i+1}:", style="Colour_Label1.TLabel")
            prio_rads_label.grid(row=i, column=0, padx=10, pady=10)

            prio_rads_entry = ttk.Entry(frame1, width=20)
            prio_rads_entry.grid(row=i, column=1, padx=10, pady=10)

            self.prio_rads_entries.append(prio_rads_entry)  # Save reference

    def get_prio_rads_data(self):
        """
        Retrieve data from the dynamically created structure entry fields.
        """
        data = [entry.get() for entry in self.prio_rads_entries]
        print("Entered prio_rads:", data)
        self.prio_rads_window.destroy()
        self.prio_rads_window = None
        return data


    def Insert_Priority_Str(self):
        if self.priority_str_window is None:
            self.priority_str_window = tk.Toplevel(self.root, padx=10, pady=10)
            self.priority_str_window.title("Priority Structure")
            self.priority_str_window.geometry("600x500")
            self.priority_str_window.configure(background="lightblue")

            frame = ttk.Frame(self.priority_str_window, style="Colour_Frame.TFrame")
            frame.grid(row=0, column=0, sticky=tk.W)

             # Scrollable frame for dynamic fields
            canvas = tk.Canvas(self.priority_str_window, background="lightblue", height=300, width=560)
            scrollbar = ttk.Scrollbar(self.priority_str_window, orient=tk.VERTICAL, command=canvas.yview)
            scrollable_frame = ttk.Frame(canvas, style = "Colour_Frame.TFrame")

            # Configure scroll region
            scrollable_frame.bind(
                "<Configure>",
                lambda e: canvas.configure(scrollregion=canvas.bbox("all")),
            )
            canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
            canvas.configure(yscrollcommand=scrollbar.set)

            # Place the scrollable frame and scrollbar
            canvas.grid(row=1, column=0, sticky=tk.W + tk.E)
            scrollbar.grid(row=1, column=1, sticky=tk.N + tk.S)

            frame1 = ttk.Frame(self.priority_str_window, style="Colour_Frame.TFrame")
            frame1.grid(row=2, column=0)

            label = ttk.Label(frame, text = "Number of Structures:", style = "Colour_Label.TLabel")
            label.grid(row = 0, column = 0 , padx = 10, pady = 10)
            str_entry = ttk.Entry(frame, width = 10 )
            str_entry.grid(row = 0, column = 1)

            Insert_button = ttk.Button(frame, text = "Insert", command =lambda:self.generate_priority_str_fields(scrollable_frame, str_entry))
            Insert_button.grid(row = 0, column = 2, padx = 10, pady = 10 )

            close_button = ttk.Button(frame1, text = "DONE", command = self.get_priority_str_data)
            close_button.grid(row = 0, column = 0, pady = 10, columnspan=2 )

    def generate_priority_str_fields(self, frame1, str_entry):
        """
        Generate labels and entry fields dynamically based on user input.
        """
        try:
            str_number = int(str_entry.get())  # Get the number of structures
        except ValueError:
            # Handle invalid input gracefully
            tk.messagebox.showerror("Invalid Input", "Please enter a valid number.")
            return

        # Clear previous fields
        for widget in frame1.winfo_children():
            widget.destroy()

        for i in range(str_number):
            str_label = ttk.Label(frame1, text=f"Structure {i+1}:", style="Colour_Label1.TLabel")
            str_label.grid(row=i, column=0, padx=10, pady=10)

            struc_entry = ttk.Entry(frame1, width=50)
            struc_entry.grid(row=i, column=1, padx=10, pady=10)

            self.prio_structure_entries.append(struc_entry)  # Save reference

    def get_priority_str_data(self):
        """
        Retrieve data from the dynamically created structure entry fields.
        """
        data = [entry.get() for entry in self.prio_structure_entries]
        print("Entered Structures:", data)
        self.priority_str_window.destroy()
        self.priority_str_window = None
        return data

    def get_keywds(self):
        "default values"
        checksym = 0
        nmbond = 0
        symm = 1
        main_bond = []
        method_type = self.update_method_type()
        if method_type == 'Chem inst':
            chinst = 1
        elif method_type == 'Rumer':
            chinst = 0

        cheminst = self.cheminst_str_type_read()
        if cheminst == 'Symmetry':
            symm = 1
        elif cheminst == 'Asymmetry':
            symm = 0
        elif cheminst == 'Checksymm':
            checksym = 1


        sotype = self.symmetry_set_order_type_read()
        if sotype == 'Quality-Arrange':
            set_order = 0
        elif sotype == 'Small-to-Big':
            set_order = 1
        elif sotype == 'Big-to-Small':
            set_order = 2


        settype = self.ChemInst_set_type_read()
        if settype == 'Single Set':
            nset = 0
        elif settype == 'All Best Sets':
            nset = 1
        elif settype == 'All Sets':
            nset = 2
        elif settype == 'Eq Bond':
            nset = 4

        mout = self.mout_number

        overlap = self.get_str_ovlp_keywd()
        if overlap == 'yes':
            ovlp = 1
        else:
            ovlp = 0

        itb = self.get_IAB_Priority()
        nnb = self.get_NAB_Priority()
        syb = self.get_SBB_Priority()
        mnbond = self.get_PDB_Priority()
        radical = self.get_PDR_Priority()

        if itb == 'None':
            itbp = 0
        else:
            itbp = int(itb)

        if nnb == 'None':
            nnbp = 0
        else:
            nnbp = int(nnb)

        if syb == 'None':
            sybp = 0
        else:
            sybp = int(syb)

        if mnbond == 'None':
            mnbondp = 0
        else:
            mnbondp = int(mnbond)

        if radical == 'None':
            radicalp = 0
        else:
            radicalp = int(radical)

        if mnbond != 'None':
            nmbond = self.bond_number
            main_bond = self.PDB_data
        
        return (int(chinst), int(symm), int(checksym), 
                int(set_order), int(nset), int(mout), int(ovlp), int(itbp), 
                int(nnbp), int(sybp), int(mnbondp), int(radicalp), int(nmbond))

class Run_Fort:
    def __init__(self,root, ctrl_class):
        self.root = root
        self.ctrl_class = ctrl_class
    def get_keywds(self, keywd_class):
        self.chinst, self.symm, self.checksym, self.set_order, self.nset,\
                self.mout, self.ovlp, self.itb, self.nnb, self.syb, self.mnbond, \
                self.radical, self.nmbond=keywd_class.get_keywds()
#        print('checksym_type',type(self.checksym))
#        print('i am in get keywds',type(self.chinst), type(self.symm), type(self.checksym), type(self.set_order), type(self.nset), type(self.mout),\
#                type(self.ovlp), type(self.itb), type(self.nnb), type(self.syb), type(self.mnbond), type(self.radical), type(self.nmbond))
#        symm_str.get_splkeywds(self.chinst, self.symm, self.checksym, self.set_order, self.nset, self.mout, self.ovlp, self.itb, self.nnb,\
#                self.syb, self.mnbond, self.radical, self.nmbond)

    def get_orbs(self, orb_class):
        self.atoset, self.norbsym, self.active, self.atn, self.orbsym = orb_class.get_orbital_matrices()
        print('orbs_datai,atoset',self.atoset)
        print('norbsym',self.norbsym)
        print('active',self.active)
        print('atn',self.atn)
        print('orbsym',self.orbsym)
#        print('orbs_data_shape',self.atoset.shape, self.norbsym.shape, self.active.shape, self.atn.shape, self.orbsym.shape)
        #try:
        #    self.atoset_py = np.zeros((200,20), dtype=int)    
        #    self.atoset_py[:len(atoset)]=atoset
#       #     self.atoset_row, self.atoset_col = self.atoset.shape

#       #     print('atoset_row, atoset_col',self.atoset_row, self.atoset_col)
        #    self.norbsym_py = np.zeros(50, dtype=int)  
        #    self.norbsym_py[:len(norbsym)] = norbsym  
#       #     self.norbsym_size = len(self.norbsym)

        #    self.active_py = np.zeros(30, dtype=int)    
        #    self.active_py[:len(active)] = active
#       #     self.active_size = len(self.active)

        #    self.atn_py = np.zeros(50, dtype=int)         
        #    self.atn_py[:len(atn)] = atn   
#       #     self.atn_size = len(self.atn)

        #    self.orbsym_py = np.zeros((20, 20), dtype=int)   
        #    self.orbsym_py[:len(orbsym)] = orbsym
#       #     self.orbsym_row, self.orbsym_col = self.orbsym.shape
        #except ValueError as e:
        #    raise ValueError(f"Error in input data conversion: {e}")
        #print('orbs_data',self.atoset, self.norbsym, self.active, self.atn, self.orbsym)
#        symm_str.get_orbs_info(self.atoset, self.norbsym, self.active, \
#                self.atn, self.orbsym, self.atoset_row, self.atoset_col, self.norbsym_size, \
#                self.active_size, self.atn_size, self.orbsym_row, self.orbsym_col)

    def get_geometry(self):
#        global readgeo, GlobVar.geometry_inserted
        if not GlobVar.geometry_inserted:
            messagebox.showerror("Invalid Geometry", "please insert the geometry")
            return
        symat, coordx, coordy, coordz, symatno = GlobVar.readgeo.get_geometry_data()
        # Convert each to a numpy array
        self.symat_py = np.zeros(20, dtype="U5")
        self.symat_py[:len(symat)]=symat
        self.coordx_py = np.zeros(100, dtype=np.float64)
        self.coordx_py[:len(coordx)]=coordx
        self.coordy_py = np.zeros(100, dtype=np.float64)
        self.coordy_py[:len(coordy)]=coordy
        self.coordz_py = np.zeros(100, dtype=np.float64)
        self.coordz_py[:len(coordz)]=coordz
        self.symatno_py = np.zeros(20, dtype=np.float64)
        self.symatno_py[:len(symatno)] = symatno

#        self.size_array=len(self.symat)

#        symm_str.get_geometry_info(self.symat, self.coordx, self.coordy, self.coordz, self.symatno, self.size_array)
        print('geometry:',self.symat_py, self.coordx_py, self.coordy_py, self.coordz_py, self.symatno_py)
    
    def get_ctrl_keywds(self, ctrl_keywds):
        self.geometry_unit, self.nao, self.nae, self.nmul, self.output_file_name = ctrl_keywds.get_ctrl_keywds()
        symm_str.get_ctrl_inputs(self.geometry_unit, self.nao, self.nae, self.nmul, self.output_file_name, self.chinst, self.symm,\
                self.checksym, self.set_order, self.nset, self.mout, self.ovlp, self.itb, self.nnb,\
                self.syb, self.mnbond, self.radical, self.nmbond, self.symat_py, self.coordx_py, self.coordy_py,\
                self.coordz_py, self.symatno_py, self.atoset, self.norbsym, self.active, self.atn, self.orbsym)

#        print(geometry_unit, nao, nae, nmul, output_file_name)
        
    def share_input_data(self):
        symm_str.get_input_data(
                str(self.geometry_unit), 
                int(self.nao), 
                int(self.nae), 
                int(self.nmul), 
                str(self.output_file_name),
                int(self.atoset_row), 
                int(self.atoset_col), 
                int(self.norbsym_size),
                int(self.active_size), 
                int(self.atn_size), 
                int(self.orbsym_row), 
                int(self.orbsym_col), 
                self.atoset, 
                self.norbsym, 
                self.active, 
                self.atn, 
                self.orbsym, 
                self.symat, 
                self.coordx, 
                self.coordy, 
                self.coordz,
                self.symatno, 
                int(self.size_array), 
                int(self.chinst), 
             #   int(self.symm), 
             #   int(self.checksym), 
             #   int(self.set_order), 
             #   int(self.nset),
             #   int(self.mout), 
             #   int(self.ovlp), 
             #   int(self.itb), 
             #   int(self.nnb), 
             #   int(self.syb), 
             #   int(self.mnbond), 
             #   int(self.radical), 
             #   int(self.nmbond)
            )






class class_manager:
    def __init__(self,root, inputc):
        self.root = root
        self.inputc = inputc

    def finish(self):
        sys.exit()
    
    def create_Ctrl_Input(self):
        self.inputc.create_ctrl_pans()
            
    def create_orb(self, keywd_button):
        self.inputo=Orb_Input( self.root, keywd_button)
    
    def create_keywd(self, Run_button):
#        global multiplicity, type_orb_count
        self.input_keywd = Keywd_Input(self.root, Run_button)
        self.input_keywd.create_keywd_pane()
    
    def Run_fort_subs(self):
        self.runfort = Run_Fort(self.root, self.inputc)
        self.runfort.get_keywds(self.input_keywd)
        self.runfort.get_orbs(self.inputo)
        self.runfort.get_geometry()
        self.runfort.get_ctrl_keywds(self.inputc)
#        self.runfort.share_input_data()

    

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Input File Creator")
    root.geometry("550x750")
    root.configure(background="lightblue")
    inputc = Ctrl_Input(root)
    manager = class_manager(root, inputc)

    style_colour_frame = ttk.Style()
    style_colour_frame.configure("Colour_Frame.TFrame", background="lightblue")
    frame1 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
    frame1.grid(row=5, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
    frame2 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
    frame2.grid(row=6, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

    Run_button = ttk.Button(frame1, text="RUN", command=lambda:manager.Run_fort_subs())
    Run_button.grid(row=1, column=0, padx=5, pady=10)
    Run_button.config(state=tk.DISABLED)  # Initially disable the button
    inputc.Run_button = Run_button

    keywd_button = ttk.Button(frame1, text="KEYWDS", command=lambda:manager.create_keywd( Run_button))
    keywd_button.grid(row=0, column=1, padx=5, pady=10)
#    keywd_button.config(state=tk.DISABLED)  # Initially disable the button
    inputc.keywd_button = keywd_button

    orbital_button = ttk.Button(frame1, text="ORBITALS", command=lambda:manager.create_orb( keywd_button))
    orbital_button.grid(row=0, column=0, padx=5, pady=10)
    orbital_button.config(state=tk.DISABLED)  # Initially disable the button
    inputc.orbital_button = orbital_button

    close_button = ttk.Button(frame2, text="FINISH", command=manager.finish)
    close_button.grid(row=0, column=1, pady=10)

    root.mainloop()
