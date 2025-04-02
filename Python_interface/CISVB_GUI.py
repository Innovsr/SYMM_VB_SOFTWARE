###################################################################################################
'''
 python GUI interface of Chemically Intuitive & Symmetric Valence Bond structure set generation
 code to interact with the users; This interface is used to get control keywords (Class
 CtrlInput): Chemical Formula (validate other inputs & create output file name), Number of Active 
 orbitals (nao), Number of Active Electrons (nae), Multiplicity (mult); Geometry of molecule
 (Class ReadGeo) : there are two options to insert geometry (1: can be brows to upload a file,
 2: can insert manually), The geometry can be inserted in two different units Bohr or Angstrom;
 Orbital informations (class Orb_Input): User need to specify the number of the active atoms
 (according to geometry) associated with active orbitals. User also need to specify type of the
 orbitals (sig or pi(x) or pi(y) or pi(z)), they can also provide fragments (the cluster of
 atomic orbitals such as s, px, py, pz, dxx, dxy, dyy etc. over which the orbitals delocalised);
 finally user need to specify spatial keywords (class Keywd_Input): method type (Chem_Inst:
 provide Chemically meaningfull VB structure set or Rumer: provide Rumer VB structure set),
 Rumer Set type (Single_Rumer_set: Rumer set corresponding to the present orbital order or
 All_Rumer_Set: Rumer sets corresponding to all possible unique permutations of the orbital
 order); Structure Type (Covalent or Ionic or Both); Chem Inst Set Type (Single set: One Best
 Chemically meaning ful set, All Best Sets: All highesr quality sets, All Sets: all possible
 sets with their quality scores, Eq Bond: Sets where all bonds appear same number of time );
 Chem_Inst Str Type ( Symmetry: sets will be symmetric, Asymmetry: there may not be any
 symmetry in the sets, Checksymm: Check if the set or sets are symmetric or not); Symmetric
 group arrangement (Quality-Arrange: Symmetric groups are arranged according to chemical
 quality scores, Big-to-Small: big group to small group, small-to-big: arrange small group to
 big group); Estimate Overlap (Yes: Estimate the overlap of the structures of a set or No:);
 Decide priorities of the Chemical Qualities: (The Chemical quality of a structures has been
 measured depending of 5 criteria which are, 1) IAB : Presence of inter atomic bonds in the
 structure -Negative quality, asi less as this type of bonds prefered 2) NAB: Near atomic Bonds:
 - positive quality -as many as this type of bond is prefered, 3) SBB: Symmetry Breaking Bonds:
 - negative quality, PDB: Pre-defined  bonds, user can provide their prefered bonds, taken as
 positive quality, PDR: Pre-defined radicals, for systems with radicals, user can provide there
 prefered radical positions or orbitals, taken as positive quality), User can choose the
 priorities of these qualities; Priority Structures (user can put some structures they prefered
 to have in the set, it is possible only if the provided structures are linearly independent);
 Maximum Number of Output File: if All sets are selected then it might be important for larger
 system to specify output set numbers, by default one output set can have only 75000 sets. After
 getting all the data, class Run_Fort send them to Fortran subroutines to generate set of VB
 structures and store them in the output file.
 '''
##################################################################################################

import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
import re
import os
import sys
import numpy as np
import Pmw
import symm_str
import tkinter.font as tkFont
import math
from collections import Counter


class GlobVar:
    "Class Contains Global Variables"
    num_orbital = None          # number of active localised VB orbitals taken in the calculation.
    num_electron = None         # number of active lectrons taken in the calculations.
    multiplicity = None         # multiplicity = 2*S + 1/2 , S = Total Spin of the system.
    type_orb_count = None       # counts types of orbitas present (types: sigma, px, py, pz).
    num_iao = None              # number of inactive orbitals in the system.
    geometry_inserted = False   # This indicates that the geometry has not been provided.
    readgeo = None              # Instance of ReadGeo class.
    molecule_string = None      # string of atoms present in the molecule or reaction.
    orbital_input = False       # specify if orbitals are inserted or not.
    orbital_data = []           #
    set_id = 0               # output set id to store present (on show) set number
    symm_key = 0
    at_list_bold = [            # List of Atoms according to periodic table
        'H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE', 'NA', 'MG', 'AL',
        'SL', 'P', 'S', 'CL', 'AR', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NR',
        'CU', 'ZN', 'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC',
        'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', 'CE',
        'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA',
        'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA'
    ]

    @staticmethod
    def configure_styles():
        style_colour_frame = ttk.Style()
        style_colour_frame.configure("Colour_Frame.TFrame", background="lightblue")
            
        style_colour_label = ttk.Style()
        style_colour_label.configure(
                "Colour_Label.TLabel", foreground="black", background="ivory", relief="flat",
                font=("Arial", 14)
                )
            
        style_colour_label1 = ttk.Style()
        style_colour_label1.configure(
                "Colour_Label1.TLabel", foreground="black", background="lightblue", font=("Arial", 14)
                )
            
        style_colour_label2 = ttk.Style()
        style_colour_label2.configure(
                "Colour_Label2.TLabel", foreground="black", background="deepskyblue", font=("Arial", 14)
                )

        style_colour_label3 = ttk.Style()
        style_colour_label3.configure(
                "Colour_Label3.TLabel", foreground="darkblue", background="lightblue", font=("Arial", 16)
                )
            
        style = ttk.Style()
        style.configure(
                "Custom.TRadiobutton", foreground="black", background="Lightblue", relief="raised",
                font=("Arial", 14)
                )


class Ctrl_Input:
    '''
    Class create the input panes to get Controll keywords: Chemical Formula (need to validate
    other inputs & create name of output file), Number of Active orbitals, Number of Active
    Electrons, Multiplicity
    '''
    def __init__(self, root):
        self.root = root
        self.input_text = ''
        self.insert = False
        self.unit_type_entry = False
        self.file_path = None

        Pmw.initialise(root)
        self.balloon = Pmw.Balloon(root)

        self.frames = {}

        # Create Frames
        frame_info = {
                "frame":(0, 0),
                "frame1":(1, 0),
                "frame2":(2, 0),
                "frame3":(3, 0),
                "frame4":(4, 0)
                }

        for frame_name, (row, column) in frame_info.items():
            frame = ttk.Frame(root, style="Colour_Frame.TFrame")
            frame.grid(row=row, column=column, sticky=tk.W)
            self.frames[frame_name] = frame

        #self.rect_border = ttk.Frame(self.frames["frame"], relief="solid", borderwidth=2, padding=5)
        #self.rect_border.grid(row=0, column=2)

        self.unit_type = tk.StringVar(value="None")
        # read prefered file name
        self.read_filename()
        self.create_geo_unit()

        # Input Fields
        self.entries = {}

        
        group_label_info = {
                "Top_label":("frame","Basic Informations About The Molecular System",0,16),
                "bottom_label":("frame4","Input Orbital Informations & Other Keywords",5,16),
                "geometry_label":("frame2","Geometry of the Molecular System",0,12),
                "Active_space_label":("frame3","Define Active Space of the System",0,12)
                }

        for label_name, (frame, text, row, font) in group_label_info.items():
            if frame in frame_info:
                label = ttk.Label(self.frames[frame], text=text, style="Colour_Label3.TLabel",
                                  font=("Airtel", font))
                label.grid(row=row, column=0, columnspan=3, padx=25, pady=25)

        # create Labels
        label_info ={
                "brows_geo_label":("frame2","Brows to Upload Geometry", 1, 0, 
                                   "If you have the geometry saved in any dat file you can brows\n"
                                   "that file and insert the geometry here using 'Brows' button; in\n" 
                                   "the geometry file you should have 4 or 5 columns first column: \n"
                                   "Atoms, second column: atomic numbers,third column: x coordinates,\n"
                                   "fourth column: y coordinates, fifth column: z coordinates"),
                "manual_geo_label":("frame2", "Insert Geometry Manually", 2, 0,
                                    "If you wish to insert the geometry manually please click\n"
                                    "'Geometry' button in the geometry file you should have 4 or 5 \n"
                                    "columns first column: Atoms, second column: atomic numbers,third\n" 
                                    "column: x coordinates, fourth column: y coordinates, fifth column:\n"
                                    "z coordinates")
                }

        for label_name, (frame, text, row, column, tooltip) in label_info.items():
            if frame in frame_info:
                label = ttk.Label(self.frames[frame], text=text, style="Colour_Label.TLabel")
                label.grid(row=row, column=column, sticky=tk.W, padx=15, pady=5)
                self.balloon.bind(label, tooltip)

        # create buttons
        buttons_info = {
                "brows_geo_button":("frame2", "Brows", self.read_geometry, 1, 1,
                                    "If you have the geometry saved in any dat file you can brows\n"
                                    "that file and insert the geometry here using 'Brows' button; in the \n"
                                    "geometry file you should have 4 or 5 columns first column: Atoms, \n"
                                    "second column: atomic numbers, third column: x coordinates, fourth\n" 
                                    "column: y coordinates, fifth column: z coordinates"),
                "manual_geo_button":("frame2", "Geometry", self.insert_geo_manually, 2, 1, 
                                     "If you wish to insert the geometry manually please click \n"
                                     "'Geometry' button in the geometry file you should have 4 or 5 columns\n"
                                     "first column: Atoms, second column: atomic numbers,third column: x \n"
                                     "coordinates, fourth column: y coordinates, fifth column: z coordinates"),
                "insert_button":("frame3", "Insert", self.validate_and_generate, 4, 1, 
                                 "Press Insert after providing all control keywords: Without\n"
                                 " inserting these control inputs you are not able to go to the next step")
                }

        for button_name, (frame, text, command, row, column, tooltip) in buttons_info.items():
            if frame in frame_info:
                button = ttk.Button(self.frames[frame], text=text, command=command)
                button.grid(row=row, column=column, sticky=tk.W, padx=5, pady=5)
                self.balloon.bind(button, tooltip)

        self.create_ctrl_pans()


    def create_geo_unit(self):
        '''
        Its creats the radio button to get the information about the units of geometry provided.
        there are only two options: Bohr or Angstrom.
        '''
        units = ["Bohr", "Angs"]
        label = ttk.Label(self.frames["frame2"], text="Unit of the Geometry Data", style="Colour_Label.TLabel")
        label.grid(row=3, column=0, sticky=tk.W, padx=15, pady=5)
        for i, unit in enumerate(units, start=1):
            button = ttk.Radiobutton(
                self.frames["frame2"],
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
            print('geo_unit', geo_unit)
            return (geo_unit)

    def read_filename(self):
        label1 = ttk.Label(self.frames["frame1"], text="Please enter the chemical formula", style="Colour_Label.TLabel")
        label1.grid(row=0, column=0, sticky=tk.W, padx=15, pady=5)
        label2 = ttk.Label(self.frames["frame1"], text="(e.g.: C6H6)", style="Colour_Label1.TLabel")
        label2.grid(row=0, column=2, sticky=tk.W, padx=5, pady=5)
        self.balloon.bind(label1, "Enter a sequence of element symbols followed by numbers \n"
                          "to specify the amounts of desired elements (e.g., C6H6, N2O,) \n"
                          "For reaction use 'under_score' in between two reactants (e.g., O_H2)")
        self.molecule_entry = ttk.Entry(self.frames["frame1"], width=15)
        self.molecule_entry.grid(row=0, column=1, padx=5, pady=5)

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
            messagebox.showerror("Incorrect Molecular Formula", "The molecular formula is not in Correct\n"
                                 "Format Enter a sequence of element symbols followed by numbers \n"
                                 "to specify the amounts of desired elements (e.g., C6H6, N2O,) \n"
                                 "For reaction use 'under_score' in between two reactants (e.g., O_H2)")
            return

    def read_geometry(self):
        ''' Open file dialog to select a file'''
        self.file_path = filedialog.askopenfilename(title="Select Geometry File")

        if self.file_path:  # Check if a file is selected
            GlobVar.readgeo = Read_Geo(self.file_path)  # Initialize the Read_Geo class
            GlobVar.readgeo.read_geometry()       # Call the read_geometry method

            ttk.Button(
                self.frames["frame2"],
                text="View_Geometry",
                command=GlobVar.readgeo.display_geometry
            ).grid(row=1, column=2, columnspan=2)

    def insert_geo_manually(self):
        """Allows manual insertion of geometry data via Read_Geo."""
        GlobVar.readgeo = Read_Geo(self.file_path)  # Initialize the Read_Geo class
        GlobVar.readgeo.insert_geo(self.root)  # Call insert_geo method
        if GlobVar.geometry_inserted is True:
            button = ttk.Button(
                    self.frames["frame1"], text="View Geometry", command=GlobVar.readgeo.display_geometry
                    )
            button.grid(row=2, column=2, sticky=tk.W, padx=5, pady=5)

    def create_ctrl_pans(self):
        self.keywords = [
                " Number of Active Orbitals (nao) : ", " Number of Active Electrons (nae) : ",
                " Multiplicity of Active Part (nmul) : "
                ]
        for idx, key in enumerate(self.keywords):
            label = ttk.Label(self.frames["frame3"], text=key, style="Colour_Label.TLabel")
            label.grid(row=idx + 1, column=0, sticky=tk.W, padx=15, pady=5)
            entry = ttk.Entry(self.frames["frame3"], width=20)
            entry.grid(row=idx + 1, column=1, padx=5, pady=5)
            self.entries[key] = entry

    def validate_and_generate(self):
        self.ctrl_inputs = []
        if GlobVar.geometry_inserted is False:
            messagebox.showerror("Geometry Error", "Dont forget to insert the geometry of the system")
        try:
            GlobVar.molecule_string = self.molecule_entry.get()
            if not GlobVar.molecule_string:
                messagebox.showerror("Molecule Error", "Please enter molecular Formula")
                return

            Total_Electrons = self.count_Total_Electron(GlobVar.molecule_string)
            print('Total Number of Electrons', Total_Electrons)
        except ValueError:
            messagebox.showerror("Molecule Error", "Please enter molecular Formula Properly ")
        if self.unit_type_entry is False:
            messagebox.showerror(
                    "Unit Type Missing", "Please Select a Unit between Bohr and Angs (angstrom)."
                    )
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
            nao, nae, nmul = self.ctrl_inputs
            GlobVar.num_orbital = int(nao)
            GlobVar.num_electron = int(nae)
            GlobVar.multiplicity = int(nmul)
#            print('nao, nae, nmul', GlobVar.num_orbital, GlobVar.num_electron, GlobVar.multiplicity)
            if (Total_Electrons - GlobVar.num_electron) % 2 == 0:
                GlobVar.num_iao = int((Total_Electrons - GlobVar.num_electron)/2)
#                print('num_iao',GlobVar.num_iao)
            else:
                messagebox.showerror("Active Orbital Error", "The number of Active Orbitals or \n"
                                     "Number of Active Electrons or Molecular Formula is not correct")
                return
            orbital_button.config(state=tk.NORMAL)
            return True  # Validation successful

    def get_ctrl_keywds(self):
        geometry_unit = self.update_geo_unit()
        nao, nae, nmul = self.ctrl_inputs
        return (geometry_unit, nao, nae, nmul)

##############################################################################################
# Reading geometry starts here :                                                             #
##############################################################################################


class Read_Geo:
    def __init__(self, file_path):
        self.file_path = file_path
        self.geo_frame = None
        style_colour_label = ttk.Style()
        style_colour_label.configure(
                "Colour_Label.TLabel", foreground="black", background="Ivory", relief="flat",
                font=("Arial", 14)
                )
        style_colour_frame = ttk.Style()
        style_colour_frame.configure("Colour_Frame.TFrame", background="lightblue")
        self.atoms = []      # List to store atom data
        self.symat = []      # list to store atom names
        self.symatno = []    # list to store atomic numbers
        self.coordx = []     # list to store atoms x coordinates
        self.coordy = []     # list to store atoms y coordinates
        self.coordz = []     # list to store atoms z coordinates

    def Deconvert(self, value):
        ''' Check if the string contains 'D' for scientific notation '''
        if 'D' in value:
            value = value.replace('D', 'e')
        return float(value)

    def read_geometry(self):
        '''
        Browse a file and read its lines to extract atom data.
        Handles two formats of geometry files:
        - 4 columns: atom name, x, y, z coordinates
        - 5 columns: atom name, atomic number, x, y, z coordinates
        '''
        try:
            with open(self.file_path, 'r') as file:
                lines = file.readlines()
            for line in lines:
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                columns = line.split()
                if len(columns) == 4:
                    ''' 4-column format: atom name, x, y, z'''
                    atom_name = columns[0].capitalize()  # Capitalize atom name
                    x = self.Deconvert(columns[1])
                    y = self.Deconvert(columns[2])
                    z = self.Deconvert(columns[3])
                    self.atoms.append({
                        "atom": atom_name,
                        "x": x,
                        "y": y,
                        "z": z
                    })
                elif len(columns) == 5:
                    ''' 5-column format: atom name, atomic number, x, y, z'''
                    atom_name = columns[0].capitalize()  # Capitalize atom name
                    atomic_number = self.Deconvert(columns[1])
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
                    messagebox.showerror(
                            "Unknown Format", "Geometry file format is unknown Please Check help"
                            )
                    continue
#            print('atom', self.atoms)
            self.total_atoms = len(self.atoms)
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
#            print("self.symat, self.coordx, self.coordy, self.coordz, self.symatno",self.symat, \
#            self.coordx, self.coordy, self.coordz, self.symatno)
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
        # frams
        frame1 = ttk.Frame(geo_window, style="Colour_Frame.TFrame", padding=10)
        frame1.grid(row=0, column=0)
        frame3 = ttk.Frame(geo_window, style="Colour_Frame.TFrame", padding=10)
        frame3.grid(row=2, column=0)
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

        label = ttk.Label(frame1, text="Please insert the number of atoms: ", style="Colour_Label.TLabel")
        label.grid(row=0, column=0, padx=10, pady=10, sticky=tk.W)

        atom_num_entry = ttk.Entry(frame1, width=10)
        atom_num_entry.grid(row=0, column=1, padx=10, pady=10)

        def create_pan():
            self.atom_num = int(atom_num_entry.get())

            for widget in frame2.winfo_children():
                widget.destroy()
            entry_widgets.clear()

            label1 = ttk.Label(frame2, text=" Atom ", background="lightblue")
            label1.grid(row=1, column=0, padx=10)
            label2 = ttk.Label(frame2, text=" X Coordinates ", background="lightblue")
            label2.grid(row=1, column=1, padx=10)
            label3 = ttk.Label(frame2, text=" Y Coordinates ", background="lightblue")
            label3.grid(row=1, column=2, padx=10)
            label4 = ttk.Label(frame2, text=" Z Coordinates ", background="lightblue")
            label4.grid(row=1, column=3, padx=10)

            for i in range(self.atom_num):
                Atom = ttk.Entry(frame2, width=10)
                Atom.grid(row=2+i, column=0, padx=10, pady=10)

                X_coord = ttk.Entry(frame2, width=15)
                X_coord.grid(row=2+i, column=1, padx=10, pady=10)

                Y_coord = ttk.Entry(frame2, width=15)
                Y_coord.grid(row=2+i, column=2, padx=10, pady=10)

                Z_coord = ttk.Entry(frame2, width=15)
                Z_coord.grid(row=2+i, column=3, padx=10, pady=10)

                entry_widgets.append((Atom, X_coord, Y_coord, Z_coord))

        def fetch_data():
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

            self.total_atoms = self.atom_num
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
#           print("self.symat, self.coordx, self.coordy, self.coordz, self.symatno", self.symat, self.coordx,\
#           self.coordy, self.coordz, self.symatno)
            GlobVar.geometry_inserted = True

            return self.symat, self.coordx, self.coordy, self.coordz, self.symatno

        button = ttk.Button(frame1, text="Enter", command=create_pan)
        button.grid(row=0, column=3, padx=10, pady=10, sticky=tk.W)
        button = ttk.Button(frame3, text="Insert", command=fetch_data)
        button.grid(row=0, column=3, padx=10, pady=10, sticky=tk.W)
        button = ttk.Button(frame3, text="Close", command=geo_window.destroy)
        button.grid(row=1, column=3, padx=10, pady=10, sticky=tk.W)

    def display_geometry(self):
        ''' Create a Toplevel window to display the file geometry '''
        if GlobVar.geometry_inserted is True:
            geo_file_display = tk.Toplevel()
            if (GlobVar.molecule_string):
                geo_file_display.title(f"{GlobVar.molecule_string}-geometry")
            else:
                geo_file_display.title("-Geometry-")
            geo_file_display.geometry("700x500")  # Optional: Set the size of the Toplevel window

            style = ttk.Style()
            style.configure("Treeview", font=("Helvetica", 12))
            # Define the Treeview widget with columns
            columns = ("atom", "x", "y", "z")
            tree = ttk.Treeview(geo_file_display, columns=columns, show="headings")
            #tree.pack(expand=True, fill="both")

            # Define headings
            tree.heading("atom", text="Atom", anchor="center")
            tree.heading("x", text="X", anchor="center")
            tree.heading("y", text="Y", anchor="center")
            tree.heading("z", text="Z", anchor="center")

            tree.column("atom", width=60, anchor="center")
            tree.column("x", width=120, anchor="center")
            tree.column("y", width=120, anchor="center")
            tree.column("z", width=120, anchor="center")
            # Insert rows from self.atoms

            tree.insert("", tk.END, values=(" ", " ", " ", " "))
            for atom in self.atoms:
                formatted_values = (
                        atom["atom"],
                        f"{float(atom['x']):.10f}",
                        f"{float(atom['y']):.10f}",
                        f"{float(atom['z']):.10f}"
                )
                tree.insert("", tk.END, values=formatted_values)
        tree.pack(expand=True, fill="both")


        # Add a scrollbar for better usability
        scrollbar = ttk.Scrollbar(geo_file_display, orient="vertical", command=tree.yview)
        tree.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")

    def get_geometry_data(self):
        return self.symat, self.coordx, self.coordy, self.coordz, self.symatno, self.total_atoms


###################################################################################
# orbital section starts here :
###################################################################################


class Orb_Input:
    def __init__(self, root, keywd_button):
        self.num_orbital = GlobVar.num_orbital
        self.orbital_frame = None
        self.keywd_button = keywd_button
        self.atm_entry = []
        self.typ_entry = []
        self.atoset = np.zeros((200, 50), dtype=int)
        self.norbsym_py = np.zeros(50, dtype=int)
        self.orbsym = np.zeros((20, 20), dtype=int)
        self.activeatoms = np.zeros(30, dtype=int)
        self.atn_vector = np.zeros(200, dtype=int)
        self.description_orb = ""
        self.create_orbital_section()

    def create_orbital_section(self):
        if GlobVar.orbital_input is False:
            self.orbital_frame = tk.Toplevel(root, padx=5, pady=5)
            self.orbital_frame.title(f"{GlobVar.molecule_string} orbital info")
            self.orbital_frame.geometry("560x560")
            self.orbital_frame.configure(background="lightblue")
            GlobVar.orbital_input = True

            for widget in self.orbital_frame.winfo_children():
                widget.destroy()

            top_label = ttk.Label(
                    self.orbital_frame, text=f"Number of active orbitals: {self.num_orbital}",
                    style="Colour_Label.TLabel"
                    )
            top_label.grid(row=0, column=0, columnspan=3, pady=15)

            label_info={
                    'label1':("Active Orbital",0),
                    'label2':("Atom Number",1),
                    'label3':("Orbital Type",2)
                }

            for label, (text, col) in label_info.items():
                label = ttk.Label(self.orbital_frame, text=text, style='Colour_Label1.TLabel')
                label.grid(row=1, column=col, padx=5, sticky='ew')

            # Scrollable container
            container = ttk.Frame(self.orbital_frame)
            container.grid(row=2, column=0, columnspan=3, sticky=tk.W+tk.E)

            canvas = tk.Canvas(container, background="lightblue", height=370, width=510)
            scrollbar = ttk.Scrollbar(container, orient=tk.VERTICAL, command=canvas.yview)
            canvas.configure(yscrollcommand=scrollbar.set)

            frame1 = ttk.Frame(canvas, style="Colour_Frame.TFrame")
            canvas.create_window((0, 0), window=frame1, anchor="nw")

            # Configure scrollable area
            frame1.bind(
                "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
            )
            canvas.create_window((0, 0), window=frame1, anchor="nw")
            canvas.configure(yscrollcommand=scrollbar.set)

            canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
            scrollbar.pack(side=tk.RIGHT, fill=tk.Y)


            insert_button = ttk.Button(
                    self.orbital_frame, text="Insert", command=self.validate_and_store_orbital_data
                    )
            insert_button.grid(row=3, column=0, columnspan=3, pady=10)

            close_button = ttk.Button(self.orbital_frame, text="Close", command=self.destroy_orbs)
            close_button.grid(row=4, column=0, pady=10, columnspan=3)

            for i in range(self.num_orbital):
                label4 = ttk.Label(
                        frame1, text=f"active orbital {i+1} ", style="Colour_Label.TLabel"
                        )
                label4.grid(row=i+3, column=0, padx=30, pady=10, sticky=tk.W)

                # Get atom number from orbital_data if available
                atom_number = GlobVar.orbital_data[i]["atom_number"] if i < len(GlobVar.orbital_data) else ""
                self.assoatm_entry = ttk.Entry(frame1, width=10)
                self.assoatm_entry.grid(row=i+3, column=1, padx=30, pady=10)
                self.assoatm_entry.insert(0, atom_number)
                self.atm_entry.append(self.assoatm_entry)

#           creating the pane for getting associated atom number
                orbital_type = GlobVar.orbital_data[i]["orbital_type"] \
                    if i < len(GlobVar.orbital_data) else ""
                self.assotyp_entry = ttk.Entry(frame1, width=10)
                self.assotyp_entry.grid(row=i+3, column=2, padx=30, pady=10, sticky=tk.E)
                self.assotyp_entry.insert(0, orbital_type)
                self.typ_entry.append(self.assotyp_entry)

    def destroy_orbs(self):
        self.orbital_frame.destroy()
        GlobVar.orbital_input = False

    def validate_and_store_orbital_data(self):
        """Validate orbital inputs and store them in a list."""
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
            GlobVar.orbital_data.append({
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
        for data in GlobVar.orbital_data:
            i = i + 1
            if 's' in data["orbital_type"].lower():
                sig_type = sig_type + 1
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
                messagebox.showerror("Input Error", f"The type of the orbital {i} is unknown. \n"
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
        self.norbsym_py[:len(norbsym)] = norbsym
        max_len = 20

        # Pad lists to be of equal length
        pxorbs_padded = pxorbs + [0] * (max_len - len(pxorbs))
        pyorbs_padded = pyorbs + [0] * (max_len - len(pyorbs))
        pzorbs_padded = pzorbs + [0] * (max_len - len(pzorbs))
        sorbs_padded = sorbs + [0] * (max_len - len(sorbs))
        self.orbsym[0, :] = pxorbs_padded
        self.orbsym[1, :] = pyorbs_padded
        self.orbsym[2, :] = pzorbs_padded
        self.orbsym[3, :] = sorbs_padded

#        print('norbsym', self.norbsym_py)
#        print('sig_type, px_type, py_type, pz_type', sig_type, px_type, py_type, pz_type )

        # Checking how many types are non-zero
        GlobVar.type_orb_count = sum(1 for count in [sig_type, px_type, py_type, pz_type] if count > 0)

        self.atoset = self.create_matrix()
        """atoset is matrix where each row represents an atom according to the geometry, if the atom ia
        an active atom, the corresponding column get '1' otherwise '0'. and the next columns contain the
        corresponding active orbital numbers associated with that atom."""
        self.keywd_button.config(state=tk.NORMAL)
#        print('atoset',self.atoset)
#        self.orbital_button.config(state=tk.DISABLED)  # Initially disable the button
#        GlobVar.orbital_button.config(state = tk.DISABLED)

    def create_matrix(self):
        ''' Step 1: Gather atom numbers and their orbitals '''
#        print('self.orbital_data',GlobVar.orbital_data)
        atom_to_orbitals = {}
        active = []
        self.nactiveatm = 0
        for orbital_number, entry in enumerate(GlobVar.orbital_data, start=1):
            atom_number = entry["atom_number"]
            if atom_number not in atom_to_orbitals:
                atom_to_orbitals[atom_number] = []
            atom_to_orbitals[atom_number].append(orbital_number)
#        print('atom_to_orbitals',atom_to_orbitals)

        ''' Step 2: Determine the maximum atom number '''
        matrix = np.zeros((200, 20), dtype=int)

        for atom_number, orbitals in atom_to_orbitals.items():
            for col_index, orbital in enumerate(orbitals):
                matrix[atom_number - 1, col_index] = orbital + GlobVar.num_iao
                ''' num_iao = number of inactive orbitals'''

        active = np.where(matrix[:, 0] != 0)[0]+1

        self.activeatoms[:len(active)] = active
        self.nactiveatm = len(active)
        if self.nactiveatm == 0:
            messagebox.showerror("Incorrect value", "Number of active atoms is zero, please check\n"
                                 "your inputs and re-insert the values")
            return

        print(active, self.activeatoms)

        atn_vec = [sum(1 for x in row if x != 0) for row in matrix]
        self.atn_vector[:len(atn_vec)] = atn_vec

        return matrix

    def get_orbital_matrices(self):
        atoset_matrix = self.atoset
        norbsym_vector = self.norbsym_py
        active_atoms = self.activeatoms
        active_atom_num = self.nactiveatm
        atn = self.atn_vector
        orbsym_matrix = self.orbsym
        return (atoset_matrix, norbsym_vector, active_atoms, atn, orbsym_matrix, active_atom_num)


class Keywd_Input:
    def __init__(self, root):
        self.root = root
#        self.Run_Button = Run_Button
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

    def create_keywd_pane(self):
        if self.keywd_window is None:
            self.keywd_window = tk.Toplevel(self.root, padx=10, pady=10)
            self.keywd_window.title("Spatial Keyword Inputs")
            self.keywd_window.geometry("750x850")
            self.keywd_window.configure(background="lightblue")

            '''dictionary to store the frame objects'''
            self.frames = {} 
            '''dictionary to store the label objects'''
            self.labels = {}
            '''dictionary to store the keyword objects'''
            self.button = {}

            '''dictionary to store the frame names and row values'''
            frames_info = {  
                "method_type_frame": 0,
                "rumer_set_type_frame": 1,
                "str_type_frame": 2,
                "ChemInst_set_type_frame": 3,
                "cheminst_str_type_frame": 4,
                "symmetry_set_order_frame": 5,
                "ovlp_cal_frame": 6,
                "prio_label_frame": 7,
                "diff_qualities_frame": 8,
                "PDB_PDR_Button_frame": 9,
                "priority_str_frame": 10,
                "mout_frame": 11,
                "end_frame": 12
            }

            # Create frames in a loop
            for frame_name, row in frames_info.items():
                frame = ttk.Frame(self.keywd_window, style="Colour_Frame.TFrame")
                frame.grid(row=row, column=0, sticky=tk.W)
                self.frames[frame_name] = frame
           
            labels_info = {
                "method_type_label": ("method_type_frame", "Select method type: ", 0, 
                                      "Select any one from the two methods provided.\n"
                                      "Select 'Rumer' to generate a Rumer structure set,\n"
                                      "no need to provide any spatial keywords. For all\n"
                                      "other types of sets please select Chem Inst."),
    
                "rumer_set_type_label": ("rumer_set_type_frame", "Rumer Set type: ", 0,
                                         "Select any one from the two Rumer set options. \n"
                                         "Select 'Single Rumer Set' to generate one Rumer structure set,\n"
                                         "corresponding to the given order of orbitals. Select \n"
                                         "'All Rumer Sets' to get all possible unique Rumer sets."),

                "str_type_label": ("str_type_frame", "Select structure type: ", 0,
                                   "Select any one from the three types. Select 'Covalent'\n"
                                   "to generate covalent sets only. Select 'Ionic' to generate \n"
                                   "ionic structure sets only. Select both to generate both \n"
                                   "types of structure sets."),

                "ChemInst_set_type_label": ("ChemInst_set_type_frame", "Chem Inst Set type:", 0,
                                            "Select any one from the four types of Chemical\n"
                                            "insight output sets. Select 'Single Set' if you want the best\n"
                                            "chemically insightful one set, select 'All Best Sets' if you\n"
                                            "want all possible same best quality sets if available, select\n"
                                            "'All Sets' if you want all possible sets but it can be huge\n"
                                            "(adjust maximum output file accordingly). Select 'Eq Bond' if\n"
                                            "you want only equally distributed bond sets, i.e., lowest \n"
                                            "overlapped sets but less chemically meaningful. These sets \n"
                                            "could help read off negative Coulson-Chirgwin weights."),

                "ChemInst_str_type_label": ("cheminst_str_type_frame", "Cheminst Str Type:", 0,
                                            "Select any one from the three types of calculations. Select\n"
                                            "'Symmetry' if you want to have symmetric sets. Select \n"
                                            "'Asymetric' if you want to have chemical insight sets. \n"
                                            "Select 'Checksymm' if you want to check some sets are \n"
                                            "symmetric or not"),

                "symmetry_set_order_label": ("symmetry_set_order_frame", "Symmetric group arrangement:", 0,
                                             "Select any order of arrangement of the symmetric \n"
                                             "group. It could help searching a symmetric set easyer.\n" 
                                             "'Quality' arrange symmetric groups from higher quality \n"
                                             "to lower quality. other two options 'big to small' and \n"
                                             "'small to big' arrange the groups according to their sizes"),

                "ovlp_cal_label": ("ovlp_cal_frame", "Estimate Overlap", 0,
                                   "To estimate overlap among the structures \n"
                                   "in each set click yes. The default is set to 'No'"),

                "mout_label": ("mout_frame", "Maximum Number of Output Files:", 0,
                               "For a large Active Space, the total number of sets can\n"
                               "reach millions or more. Each output file can contain\n"
                               "up to 75,000 sets. By default, the number of output files\n"
                               "is set to one. If additional output files are required, please\n"
                               "specify the desired number in the entry box and press 'Enter'."),
                               
                "prio_label": ("prio_label_frame", "Decide Priorities of the Chemical Qualities:", 0,
                               "The sequence of the default priority is IAB > NAB > SBB and \n"
                               "PBU & PRU are not taken in the calculation; To change this default \n"
                               "priorities please click the button and selevct the numbers for each "),

                "IAB_label": ("diff_qualities_frame", "IAB  PRIORITY: ", 0,
                              "Intra Atomic Bond Priority. Default is set to '1'"),

                "NAB_label": ("diff_qualities_frame", "NAB  PRIORITY: ", 1,
                              "Near Atomic Bond Priority. Default is set to '2'"),

                "SBB_label": ("diff_qualities_frame", "SBB  PRIORITY: ", 2,
                              "Symmetry Breaking Bond Priority. Default is set to '3'"),

                "PDB_label": ("diff_qualities_frame", "PDB  PRIORITY: ", 3,
                              "Pre-defined Bond Priority. Default is set to \n"
                              "'None': Not taken into the calculation"),

                "PDR_label": ("diff_qualities_frame", "PDR  PRIORITY: ", 4,
                              "Pre-defined Radical Priority. Default is set to \n"
                              "'None': Not taken into the calculation")
                }
                
            for label_name, (frame_name, text, row, tooltip) in labels_info.items():
                styles = "Colour_Label.TLabel"
                if label_name in ("IAB_label", "NAB_label", "SBB_label", "PDB_label", "PDR_label"):
                    styles = "Colour_Label2.TLabel"

                if frame_name in self.frames:
                    label = ttk.Label(self.frames[frame_name], text=text, style=styles)
                    label.grid(row=row, column=0, sticky=tk.W, padx=10, pady=10)
                    self.balloon.bind(label, tooltip)
                    self.labels[label_name] = label

            self.radio_buttons_info = {
                    "Rumer_Set_Type":["Single Rumer Set", "All Rumer Sets"],
                    "tip_method":["Chem inst", "Rumer"],
                    "str_type":["Covalent", "Ionic", "Both"],
                    "cheminstsets_type":["Single Set", "All Best Sets", "All Sets", "Eq Bond"],
                    "Cheminst_str_type":["Symmetry", "Asymmetry", "Checksymm"],
                    "set_order":["Quality-Arrange", "Big-to-Small", "Small-to-Big"],
                    "ovlp_types":["Yes", "No"],
                    "IAB_option":["1", "2", "3", "4", "5", "None"],
                    "NAB_option":["1", "2", "3", "4", "5", "None"],
                    "SBB_option":["1", "2", "3", "4", "5", "None"],
                    "PDB_option":["1", "2", "3", "4", "5", "None"],
                    "PDR_option":["1", "2", "3", "4", "5", "None"]
                    }

            self.radio_buttons_variable = {
                    "Rumer_Set_Type":["rumer_set_type_frame", self.rumer_set_type, 
                                      self.Update_Rumer_Set_Type, 0],
                    "tip_method":["method_type_frame", self.method_type, self.update_method_type, 0],
                    "str_type":["str_type_frame", self.str_type, self.Update_Str_Type, 0],
                    "cheminstsets_type":["ChemInst_set_type_frame", self.ChemInst_set_type, 
                                         self.ChemInst_set_type_read, 0],
                    "Cheminst_str_type":["cheminst_str_type_frame", self.cheminst_str_type, 
                                         self.cheminst_str_type_read, 0],
                    "set_order":["symmetry_set_order_frame", self.symmetry_set_order_type, 
                                 self.symmetry_set_order_type_read, 0],
                    "ovlp_types":["ovlp_cal_frame", self.ovlp_cal, self.get_str_ovlp_keywd, 0],
                    "IAB_option":["diff_qualities_frame", self.IAB_type, self.get_IAB_Priority, 0],
                    "NAB_option":["diff_qualities_frame", self.NAB_type, self.get_NAB_Priority, 1],
                    "SBB_option":["diff_qualities_frame", self.SBB_type, self.get_SBB_Priority, 2],
                    "PDB_option":["diff_qualities_frame", self.PDB_type, self.get_PDB_Priority, 3],
                    "PDR_option":["diff_qualities_frame", self.PDR_type, self.get_PDR_Priority, 4]
                    }

            for button, (frame, variable, callback, row) in self.radio_buttons_variable.items():
                if button in self.radio_buttons_info:
                    self.create_radiobuttons(
                            self.frames[frame], self.radio_buttons_info[button], variable, callback, row
                            )

            button_info = {
                    "PDB_button":["PDB_PDR_Button_frame", " Insert PDB ", self.Insert_PDB, 0, 0,
                                  "Click the button to insert your preffered bonds"],
                    "PDR_button":["PDB_PDR_Button_frame", " Insert PDR ", self.Insert_PDR, 0, 1,
                                  "Click the button to insert your preffered bonds"],
                    "prio_struc_button":["priority_str_frame", "Insert Priorities Structures", 
                                         self.Insert_Priority_Str, 0, 0, 
                                         "If you want some structures must be present in the set \n"
                                         "you can insert them by clicking this Button;\n"
                                         "The provided structures must be linearly independent \n"
                                         "otherwise they will not be present in the set together"],
                    "insert_button":["end_frame", "Insert", self.get_keywds, 0, 0,
                                     "Click the button to insert all values"],
                    "close_button":["end_frame", "Close", self.keywd_window.destroy, 1, 0,
                                     "Click the button to closr this Keyword pane"]
                    }
            for button_name, (frame, text, command, row, column, tooltip) in button_info.items():
                button = ttk.Button(self.frames[frame], text=text, command=command)
                button.grid(row=row, column=column, padx=10, pady=10, sticky=tk.W)
                self.balloon.bind(button, tooltip)
                self.button[button_name] = button

                if button_name in ("PDB_button", "PDR_button"):
                    button.config(state=tk.DISABLED)  # Initially disable the button


    def create_radiobuttons(self, frame, options, variable, callback, row):
        ''' Create radio buttons; check conditions and disabled specific buttons accordingly '''
        for i, option in enumerate(options, start=1):
            button = ttk.Radiobutton(
                frame,
                text=option,
                value=option,
                variable=variable,
                command=callback,
                style="Custom.TRadiobutton"
            )
            button.grid(row=row, column=i, padx=10, pady=10)
            self.disable_widgets_in_frame(self.frames["rumer_set_type_frame"])
            if GlobVar.multiplicity == 1:
                if options is self.radio_buttons_info["PDR_option"]:
                    button.config(state=tk.DISABLED)
            if self.type_orb_count == 1:
                if options is self.radio_buttons_info["SBB_option"]:
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
            self.button["PDB_button"].config(state=tk.NORMAL)
        return (PDB_Priority)

    def get_PDR_Priority(self):
        PDR_Priority = self.PDR_type.get()
        print('PDR_Priority', PDR_Priority)
        if PDR_Priority is not None and self.multiplicity != 1:
            self.button["PDR_button"].config(state=tk.NORMAL)
        return (PDR_Priority)

    def update_method_type(self):
        method_type = self.method_type.get()
        if method_type:
            self.method_type_entry = True
            if method_type == 'Rumer':
                self.enable_widgets_in_frame(self.frames["rumer_set_type_frame"])
                self.disable_widgets_in_frame(self.frames["ChemInst_set_type_frame"])
                self.disable_widgets_in_frame(self.frames["cheminst_str_type_frame"])
                self.disable_widgets_in_frame(self.frames["symmetry_set_order_frame"])
                self.disable_widgets_in_frame(self.frames["ovlp_cal_frame"])
                self.disable_widgets_in_frame(self.frames["prio_label_frame"])
                self.disable_widgets_in_frame(self.frames["diff_qualities_frame"])
                self.disable_widgets_in_frame(self.frames["PDB_PDR_Button_frame"])
                self.disable_widgets_in_frame(self.frames["priority_str_frame"])
                self.disable_widgets_in_frame(self.frames["mout_frame"])
            elif method_type == 'Chem inst':
                self.disable_widgets_in_frame(self.frames["rumer_set_type_frame"])
                self.enable_widgets_in_frame(self.frames["ChemInst_set_type_frame"])
                self.enable_widgets_in_frame(self.frames["cheminst_str_type_frame"])
                self.enable_widgets_in_frame(self.frames["symmetry_set_order_frame"])
                self.enable_widgets_in_frame(self.frames["ovlp_cal_frame"])
                self.enable_widgets_in_frame(self.frames["prio_label_frame"])
                self.enable_widgets_in_frame(self.frames["diff_qualities_frame"])
                self.enable_widgets_in_frame(self.frames["PDB_PDR_Button_frame"])
                self.enable_widgets_in_frame(self.frames["priority_str_frame"])
                self.enable_widgets_in_frame(self.frames["mout_frame"])
#            print('method_type',method_type)
            return (method_type)

    def disable_widgets_in_frame(self, frame):
        """
        Recursively disable all widgets and frames inside the given frame.
        """
        for widget in frame.winfo_children():
            widget.configure(state="disabled")

    def enable_widgets_in_frame(self, frame):
        """
        Recursively enable all widgets and frames inside the given frame.
        """
        for widget in frame.winfo_children():
            widget.configure(state="normal")

    def Update_Rumer_Set_Type(self):
        rumer_set_num = self.rumer_set_type.get()
        if rumer_set_num:
            self.rumer_set_num_entry = True
            print('rumer_set_num', rumer_set_num)
            return (rumer_set_num)

    def Update_Str_Type(self):
        structure_type = self.str_type.get()
        if structure_type:
            self.structure_type_entry = True
            print('structure_ype', structure_type)
            return (structure_type)

    def ChemInst_set_type_read(self):
        set_type = self.ChemInst_set_type.get()
        if set_type:
            self.ChemInst_set_type_entry = True
            print('set_type', set_type)
            return (set_type)

    def get_maximum_num_output(self, mout_entry):
        try:
            self.mout_number = int(mout_entry.get())
            return (self.mout_number)
        except ValueError:
            # Handle invalid input
            tk.messagebox.showerror("Invalid Input", "Please enter a valid number.")
            return
#        print('mout_number', self.mout_number)

    def get_str_ovlp_keywd(self):
        ovlp_calc = self.ovlp_cal.get()
#        print('ovlp_cal', ovlp_calc)
        return (ovlp_calc)

    def cheminst_str_type_read(self):
        cheminst_type = self.cheminst_str_type.get()
        if cheminst_type:
            self.cheminst_str_type_entry = True
#            print('cheminst_type',cheminst_type)
        return (cheminst_type)

    def symmetry_set_order_type_read(self):
        set_order_type = self.symmetry_set_order_type.get()
        if set_order_type:
            self.symmetry_set_order_type_entry = True
#            print('set_order_type',set_order_type)
            return (set_order_type)

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
            scrollable_frame = ttk.Frame(canvas, style="Colour_Frame.TFrame")

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

            label = ttk.Label(frame, text="Number of Structures:", style="Colour_Label.TLabel")
            label.grid(row=0, column=0, padx=10, pady=10)
            prio_bond_entry = ttk.Entry(frame, width=10)
            prio_bond_entry.grid(row=0, column=1)

            Insert_button = ttk.Button(
                    frame, text="Insert", command=lambda: self.generate_prio_bond_fields
                    (scrollable_frame, prio_bond_entry)
                    )
            Insert_button.grid(row=0, column=2, padx=10, pady=10)

            close_button = ttk.Button(frame1, text="DONE", command=self.get_prio_bond_data)
            close_button.grid(row=0, column=0, pady=10, columnspan=2)

    def generate_prio_bond_fields(self, frame1, prio_bond_entry):
        """
        Generate labels and entry fields dynamically based on user input.
        """
        try:
            self.bond_number = int(prio_bond_entry.get())
        except ValueError:
            tk.messagebox.showerror("Invalid Input", "Please enter a valid number.")
            return

        # Clear previous fields
        for widget in frame1.winfo_children():
            widget.destroy()

        for i in range(self.bond_number):
            prio_bond_label = ttk.Label(
                    frame1, text=f"Pre-defined Bond {i+1}:", style="Colour_Label1.TLabel"
                    )
            prio_bond_label.grid(row=i, column=0, padx=10, pady=10)

            prio_bond_entry = ttk.Entry(frame1, width=20)
            prio_bond_entry.grid(row=i, column=1, padx=10, pady=10)

            self.prio_bond_entries.append(prio_bond_entry)

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
            scrollable_frame = ttk.Frame(canvas, style="Colour_Frame.TFrame")

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

            label = ttk.Label(frame, text="Number of radicals:", style="Colour_Label.TLabel")
            label.grid(row=0, column=0, padx=10, pady=10)
            prio_rads_entry = ttk.Entry(frame, width=10)
            prio_rads_entry.grid(row=0, column=1)

            Insert_button = ttk.Button(
                    frame, text="Insert",
                    command=lambda: self.generate_prio_rads_fields(scrollable_frame, prio_rads_entry)
                    )
            Insert_button.grid(row=0, column=2, padx=10, pady=10)

            close_button = ttk.Button(frame1, text="DONE", command=self.get_prio_rads_data)
            close_button.grid(row=0, column=0, pady=10, columnspan=2)

    def generate_prio_rads_fields(self, frame1, prio_rads_entry):
        """
        Generate labels and entry fields dynamically based on user input.
        """
        try:
            rads_number = int(prio_rads_entry.get())
        except ValueError:
            tk.messagebox.showerror("Invalid Input", "Please enter a valid number.")
            return

        # Clear previous fields
        for widget in frame1.winfo_children():
            widget.destroy()

        for i in range(rads_number):
            prio_rads_label = ttk.Label(
                    frame1, text=f"Pre-defined rad {i+1}:", style="Colour_Label1.TLabel"
                    )
            prio_rads_label.grid(row=i, column=0, padx=10, pady=10)

            prio_rads_entry = ttk.Entry(frame1, width=20)
            prio_rads_entry.grid(row=i, column=1, padx=10, pady=10)

            self.prio_rads_entries.append(prio_rads_entry)

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
            scrollable_frame = ttk.Frame(canvas, style="Colour_Frame.TFrame")

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

            label = ttk.Label(frame, text="Number of Structures:", style="Colour_Label.TLabel")
            label.grid(row=0, column=0, padx=10, pady=10)
            str_entry = ttk.Entry(frame, width=10)
            str_entry.grid(row=0, column=1)

            Insert_button = ttk.Button(
                    frame, text="Insert",
                    command=lambda: self.generate_priority_str_fields(scrollable_frame, str_entry)
                    )
            Insert_button.grid(row=0, column=2, padx=10, pady=10)

            close_button = ttk.Button(frame1, text="DONE", command=self.get_priority_str_data)
            close_button.grid(row=0, column=0, pady=10, columnspan=2)

    def generate_priority_str_fields(self, frame1, str_entry):
        """
        Generate labels and entry fields dynamically based on user input.
        """
        try:
            str_number = int(str_entry.get())
        except ValueError:
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
        """Retrieve and process default keyword values."""
        def map_string_to_int(value, mapping, default=0):
            """Helper function to map string values to integers."""
            return mapping.get(value, default)

        checksym, nmbond, symm = 0, 0, 1

        method_mapping = {'Chem inst': 1, 'Rumer': 0}
        chinst = map_string_to_int(self.update_method_type(), method_mapping)

        cheminst_mapping = {'Symmetry': 1, 'Asymmetry': 0, 'Checksymm': 1}
        cheminst = self.cheminst_str_type_read()
        symm = cheminst_mapping.get(cheminst, symm)
        checksym = 1 if cheminst == 'Checksymm' else 0

        set_order_mapping = {'Quality-Arrange': 0, 'Small-to-Big': 1, 'Big-to-Small': 2}
        set_order = map_string_to_int(self.symmetry_set_order_type_read(), set_order_mapping)

        str_type_mapping = {'Both': 1, 'Covalent': 2, 'Ionic': 3}
        strtype = map_string_to_int(self.Update_Str_Type(), str_type_mapping)

        settype_mapping = {'Single Set': 0, 'All Best Sets': 1, 'All Sets': 2, 'Eq Bond': 4}
        nset = map_string_to_int(self.ChemInst_set_type_read(), settype_mapping)

        mout = self.mout_number

        overlap = 1 if self.get_str_ovlp_keywd() == 'yes' else 0

        def get_priority(value):
            """Convert priority values from string to integer."""
            return 0 if value == 'None' else int(value)

        itbp = get_priority(self.get_IAB_Priority())
        nnbp = get_priority(self.get_NAB_Priority())
        sybp = get_priority(self.get_SBB_Priority())
        mnbondp = get_priority(self.get_PDB_Priority())
        radicalp = get_priority(self.get_PDR_Priority())

        if self.get_PDB_Priority() != 'None':
            nmbond = self.bond_number

        Run_button.config(state=tk.NORMAL)

        return (chinst, symm, checksym, strtype, set_order, nset, mout, overlap,
                itbp, nnbp, sybp, mnbondp, radicalp, nmbond)


class Run_Fort:
    def __init__(self, root, ctrl_class):
        self.root = root
        self.ctrl_class = ctrl_class

    def get_keywds(self, keywd_class):
        self.chinst, self.symm, self.checksym, self.strtype, self.set_order, self.nset, \
                self.mout, self.ovlp, self.itb, self.nnb, self.syb, self.mnbond, \
                self.radical, self.nmbond = keywd_class.get_keywds()
        GlobVar.symm_key = self.symm

    def get_orbs(self, orb_class):
        self.atoset, self.norbsym, self.active, self.atn, self.orbsym, self.actv_atm_num\
                = orb_class.get_orbital_matrices()
        print('orbs_datai,atoset', self.atoset)
        print('norbsym', self.norbsym)
        print('active', self.active)
        print('atn', self.atn)
        print('orbsym', self.orbsym)
        print('actv_atm_num', self.actv_atm_num)

    def get_geometry(self):
        if not GlobVar.geometry_inserted:
            messagebox.showerror("Invalid Geometry", "please insert the geometry")
            return
        symat, coordx, coordy, coordz, symatno, self.totalatoms = GlobVar.readgeo.get_geometry_data()
        # Convert each to a numpy array
        self.symat_py = np.zeros(100, dtype="U5")
        self.symat_py[:len(symat)] = symat
        self.coordx_py = np.zeros(100, dtype=np.float64)
        self.coordx_py[:len(coordx)] = coordx
        self.coordy_py = np.zeros(100, dtype=np.float64)
        self.coordy_py[:len(coordy)] = coordy
        self.coordz_py = np.zeros(100, dtype=np.float64)
        self.coordz_py[:len(coordz)] = coordz
        self.symatno_py = np.zeros(100, dtype=np.float64)
        self.symatno_py[:len(symatno)] = symatno
#        print('geometry:',self.symat_py, self.coordx_py, self.coordy_py, self.coordz_py, self.symatno_py)

    def get_ctrl_keywds(self, ctrl_keywds):
        self.geometry_unit, self.nao, self.nae, self.nmul = ctrl_keywds.get_ctrl_keywds()

    def share_input_data(self, output_folder):
        print(self.geometry_unit)
        print(self.nao)
        print(self.nae)
        print(self.nmul)
        print(GlobVar.molecule_string)
        print(self.chinst, self.symm)
        print(self.checksym)
        print(self.set_order)
        print(self.nset)
        print(self.mout)
        print(self.ovlp)
        print(self.itb)
        print(self.nnb)
        print(self.syb)
        print(self.mnbond)
        print(self.radical)
        print(self.nmbond)
        print(self.symat_py)
        print(self.coordx_py)
        print(self.coordy_py)
        print(self.coordz_py)
        print(self.symatno_py)
        print(self.atoset)
        print(self.norbsym)
        print(self.active)
        print(self.atn)
        print(self.orbsym)
        print(self.strtype)
        print(self.totalatoms)
        print(GlobVar.num_iao)
        print(self.actv_atm_num)
        print(output_folder)
        symm_str.get_ctrl_inputs(
                self.geometry_unit,
                self.nao,
                self.nae,
                self.nmul,
                GlobVar.molecule_string,
                self.chinst, 
                self.symm,
#                self.checksym,
                self.set_order,
                self.nset,
                self.mout,
                self.ovlp,
                self.itb,
                self.nnb,
                self.syb,
                self.mnbond,
                self.radical,
                self.nmbond,
                self.symat_py,
                self.coordx_py,
                self.coordy_py,
                self.coordz_py,
                self.symatno_py,
                self.atoset,
                self.norbsym,
                self.active,
                self.atn,
                self.orbsym,
                self.strtype,
                self.totalatoms,
                GlobVar.num_iao,
                self.actv_atm_num,
                output_folder
                )
        return (self.mout)


class Output:
    def __init__(self, root):
        self.root = root

    def load_structure_file(self, fname):
        Output_window = tk.Toplevel(self.root)
        Output_window.title("All Structures")
        Output_window.geometry("1050x950")
        Output_window.configure(background="lightblue")
        structures = []
        various_qualities = []
        overall_qualities = []
        rumers = []
        sls = []

        self.frames = {}

        frame_info = {
                "frame":("lightblue", 0, 0),
                "text_frame":("white",1, 0),
                "button_frame":("lightblue",2,0),
                "set_frame":("white",3,0),
                "info_frame":("white",3,1),
                "info_button_frame":("lightblue",4,0)
                }
        for frame_name,(bg, row, column) in frame_info.items():
            frame = tk.Frame(Output_window, bg=bg)
            frame.grid(row=row, column=column, padx=10, pady=10)
            self.frames[frame_name] = frame

#        self.frames["info_button_frame"].config(style="Colour_Frame.TFrame")
        scrollbar = tk.Scrollbar(self.frames["text_frame"])
#        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        scrollbar.grid(row=0, column=1, sticky="ns")

        text_widget = tk.Text(
            self.frames["text_frame"], wrap=tk.WORD, bg="white", fg="black",
            width=120, height=20, yscrollcommand=scrollbar.set
        )
        text_widget.grid(row=0, column=0, padx=5, pady=5, sticky="w")
#        text_widget.pack(expand=True, fill=tk.BOTH)
        scrollbar.config(command=text_widget.yview)
        font = tkFont.Font(family="Helvetica", size=14)

        nlp = GlobVar.num_orbital - GlobVar.num_electron
        nao = GlobVar.num_orbital
        nmult = GlobVar.multiplicity
        sets, tot_perm_str, all_str = self.wigner(nlp, nao, nmult)

        filename = fname + '/' + 'structures.dat'

        if not os.path.exists(filename):
            text_widget.insert(tk.END, "No Structures are available\n")
            return
        else:
            i = 0
            with open(filename, "r") as file:
                lines = file.readlines()
                for line in lines:
                    if line.startswith("structure"):
                        # print('line', line)
                        cols = line.strip().split()
                        sl = cols[1]
                        sls.append(int(sl))
                        various_quality = " ".join(cols[3:8])
                        various_qualities.append(various_quality)
                        overall_quality = cols[10]
                        overall_qualities.append(int(overall_quality))
                        rumer = cols[12]
                        rumers.append(rumer)
                        if rumer == 'Rumer':
                            rumer = 'R'
                        structure = " ".join(cols[13:])
                        structures.append(structure)
                        sw = len(structure) + 10
                        i += 1
                        result = (
                                f"{i:<16}"
                                f"{structure:^{sw}}"
                                f"{various_quality:^30}"
                                f"{overall_quality:^22}"
                                f"{rumer:^20}\n\n"
                                )

                        if i == 1:
                            header1 = (
                            f"{'Sl. No.':<16}"
                            f"{'Structure':^{sw}}"
                            f"{'Various Qualities':^30}"
                            f"{'Overall Quality':^22}"
                            f"{'Rumer':^20}\n\n"
                            )

                            header2 = (
                                    f"{'-------':<16}"
                                    f"{'---------':^{sw}}"
                                    f"{'IAB NAB SBB PDB PDR':^30}"
                                    f"{'---------------':^22}"
                                    f"{'-----':^20}\n\n"
                                    )

                            text_widget.insert(tk.END, header1)
                            text_widget.insert(tk.END, header2)
                        text_widget.insert(tk.END, result)
                label_tot_cov_str = ttk.Label(
                        self.frames["frame"], text=f'Total number of covalent structure of thesystem = {all_str}', 
                        style="Colour_Label.TLabel"
                        )
                label_tot_cov_str.grid(row=0, column=0, padx=5, pady=5)
                label_tot_all_cov_str = ttk.Label(
                        self.frames["frame"], text=f'Total number of allowed covalent structure={tot_perm_str}', 
                        style="Colour_Label.TLabel"
                        )
                label_tot_all_cov_str.grid(row=0, column=1, padx=5, pady=5)
        text_widget.config(state=tk.DISABLED)

        self.tempfname = fname + '/' + 'out.temp'

#        Output_window.grid_columnconfigure(0, weight=1)
#        Output_window.grid_columnconfigure(1, weight=1)

        set_text_wt = tk.Text(self.frames["set_frame"], wrap=tk.WORD, bg="white", fg="black", width=50, height=15)
        set_text_wt.grid(row=0, column=0, padx=5, pady=5, sticky="w")

        v_set_scrollbar = tk.Scrollbar(self.frames["set_frame"], orient=tk.VERTICAL, command=set_text_wt.yview)
        v_set_scrollbar.grid(row=0, column=1, sticky="ns")

        font = tkFont.Font(family="Helvetica", size=16)
        set_text_wt.configure(yscrollcommand=v_set_scrollbar.set, font=font)

        self.set_info_wt = tk.Text(self.frames["set_frame"], wrap=tk.WORD, bg="white", fg="black", width=30, height=15)
        self.set_info_wt.grid(row=0, column=2, padx=5, pady=5, sticky="w")

        v_info_scrollbar = tk.Scrollbar(self.frames["set_frame"], orient=tk.VERTICAL, command=self.set_info_wt.yview)
        v_info_scrollbar.grid(row=0, column=3, sticky="ns")

        font = tkFont.Font(family="Helvetica", size=16)
        self.set_info_wt.configure(yscrollcommand=v_info_scrollbar.set, font=font)

#        self.frames["set_frame"].grid(row=3, column=0, sticky="w")
#        self.frames["info_frame"].grid(row=3, column=1, sticky="ew")

        print('out.temp file path', self.tempfname)
        if not os.path.exists(self.tempfname):
            messagebox.showerror(
                    'Error', 'No Output File has been found'
                    )
            #set_text_wt.insert(tk.END, "No Structures are available\n")
            return
        else:
            with open(self.tempfname, "r") as file:
                try:
                    lines = file.readlines()
                    print('lines',lines)
                except: 
                    messagebox.showerror(
                            'Error', 'No independent set has been found, try after few munits'
                            )
                    return

        flag1 = 1
        flag2 = 2
        flag3 = 3
        view_set_button = tk.Button(
                self.frames["button_frame"], 
                text='View Set', 
                command=lambda: self.view_set(flag1, structures, various_qualities, overall_qualities,
                                              rumers, sls, set_text_wt, lines)
                )
        view_set_button.grid(row=0, column=2, padx=10, pady=10)

        view_nextset_button = tk.Button(
                self.frames["button_frame"], text='Next Set', 
                command=lambda: self.view_set(flag2, structures, various_qualities, overall_qualities,
                                              rumers, sls, set_text_wt, lines)
                )
        view_nextset_button.grid(row=0, column=3, padx=10, pady=10)

        view_prevtset_button = tk.Button(
                self.frames["button_frame"], text='Prev Set', 
                command=lambda: self.view_set(flag3, structures, various_qualities, overall_qualities,
                                              rumers, sls, set_text_wt, lines)
                )
        view_prevtset_button.grid(row=0, column=1, padx=10, pady=10)


    def view_set(self, flag, structure, various_qualities, overall_qualities, rumers, sls, set_text_wt, lines):
        print('inside view_set',flag)
        if flag == 1:
            GlobVar.set_id = 1
            self.set_info_wt.delete("1.0", "end")
        elif flag == 2:
            GlobVar.set_id += 1
            self.set_info_wt.delete("1.0", "end")
        elif flag == 3:
            GlobVar.set_id -= 1
            self.set_info_wt.delete("1.0", "end")
        print('flag',GlobVar.set_id)

        set_text_wt.tag_configure("right", justify="right")
        #print('out.temp file path', self.tempfname)
        #if not os.path.exists(self.tempfname):
        #    set_text_wt.insert(tk.END, "No Structures are available\n")
        #    return
        #else:
        #    with open(self.tempfname, "r") as file:
        #        try:
        #            lines = file.readlines()
        #        except: 
        #            messagebox.showerror(
        #                    'Error', 'No independent set has been found, try after few munits'
        #                    )
        #            return
        print('lines',lines)
        for line in lines:
            if line.startswith('Set_number'):
                linear_indset = line.strip().split()

                if linear_indset[1] == str(GlobVar.set_id):
                    indset = []
                    set_text_wt.delete("1.0", "end")
                    set_text_wt.insert(
                            tk.END, "\n           Independent Set of structures  \n\n"
                            )
                    set_text_wt.insert(
                            tk.END, "\n           Set number::  " + f"{linear_indset[1]}\n\n"
                            )
                    for index in linear_indset[2:]:
                        indset.append(int(index))
                        set_text_wt.insert(
                                tk.END, f"   {sls[int(index) - 1]}:        " + structure[int(index) - 1] + "\n"
                                )
                    info_button = tk.Button(self.frames["info_button_frame"], text="Set Info", 
                           command=lambda: self.set_info(
                                various_qualities, overall_qualities, rumers, indset
                                )
                            )
                    info_button.grid(row=0, column=0, sticky=tk.W, padx=10)
                    if GlobVar.symm_key == 1:
                        symm_grp_button = tk.Button(self.frames["info_button_frame"], text="Symmetry groups",
                                                    command=lambda:self.show_symm_groups(
                                                        overall_qualities, structure
                                                        )
                                                    )
                        symm_grp_button.grid(row=0, column=1, sticky=tk.W, padx=10)
                    print('linear_indset', linear_indset)

    def show_symm_groups(self, overall_qualities, structure):
        strpergrp = Counter(overall_qualities)
        strpergrp = dict(sorted(strpergrp.items()))
        print('strpergrp', strpergrp)

    def set_info(self, various_qualities, overall_qualities, rumer, indset):
        oqt = sum(overall_qualities[i - 1] for i in indset)
        self.set_info_wt.delete("1.0", "end")
        self.set_info_wt.insert(tk.END, f"\n Overall Quality: {oqt}\n\n")
        self.set_info_wt.insert(tk.END, "IAB  NAB  SBB  PDB  PDR\n")
        self.set_info_wt.insert(tk.END, f"{'-----':<8}" 
                                f"{'-----':<8}" 
                                f"{'-----':<8}" 
                                f"{'-----':<8}" 
                                f"{'-----':<8}\n")
        for i in indset:
            a, b, c, d, e = various_qualities[i-1].strip().split()
            self.set_info_wt.insert(tk.END, f"{a:<8}" 
                                    f"{b:<8}" 
                                    f"{c:<8}" 
                                    f"{d:<8}" 
                                    f"{e:<8}""\n")
        j = 0
        for i in indset:
            if rumer[i-1] == 'R':
                j += 1

        if j == len(indset):
            self.set_info_wt.insert(tk.END, "\n The set is a Rumer Set\n")
        else:
            self.set_info_wt.insert(tk.END, "\n The set is a Chemical Insight Set \n")

        return #oqt

    def wigner(self, nlp, nao, nmult):
        ''' none = number of one electron orbital
            nlp = number of lone paires
            nao = number of active electrons
            nmult = multiplicity
        '''
        wigner = 0
        sets = 1

        ''' calculate permissible number of structre with Wigner's theorem'''
        noeo = nao - nlp
        spin = (nmult - 1)/2
        neum = nmult*math.factorial(noeo)
        denom1 = math.factorial(int((noeo/2) + spin+1))
        denom2 = math.factorial(int((noeo/2) - spin))
        wigner = int(neum/(denom1 * denom2))

        ''' Number of sets depending on lone pairs'''
        if nlp != 0:
            sets = math.comb(nao, nlp)

        ''' Calculate total number of structures'''
        term1 = math.comb(nao, nao - noeo)
        term2 = math.comb(noeo, int(2*spin))

        product_term = np.prod([(math.comb(noeo - int(2*spin) - 2*i, 2))
                                for i in range(int((noeo / 2)-spin))])
        denom = math.factorial(int((noeo / 2) - spin))

        totstr = int((term1 * term2 * product_term) / denom)

        return (sets, wigner, totstr)


class class_manager:
    def __init__(self, root, inputc):
        self.root = root
        self.inputc = inputc

    def finish(self):
        sys.exit()

    def create_Ctrl_Input(self):
        self.inputc.create_ctrl_pans()

    def create_orb(self, keywd_button):
        self.inputo = Orb_Input(self.root, keywd_button)

    def create_keywd(self, Run_button):
        self.input_keywd = Keywd_Input(self.root)
        self.input_keywd.create_keywd_pane()
        self.Run_button = Run_button

    def Run_fort_subs(self):
        self.runfort = Run_Fort(self.root, self.inputc)
        self.runfort.get_keywds(self.input_keywd)
        self.runfort.get_orbs(self.inputo)
        self.runfort.get_geometry()
        self.runfort.get_ctrl_keywds(self.inputc)

        cwd = os.getcwd()
        last_part = '_output'
        end_name = f"{GlobVar.molecule_string}{last_part}"
        self.output_folder = os.path.join(cwd, end_name)
#        print(f"Current Working Directory:{self.output_folder}")
        try:
            os.makedirs(self.output_folder, exist_ok=True)
            '''`exist_ok=True` avoids error if the folder already exists'''
        except OSError as e:
            tk.messagebox.showerror("Error", f"{e}")
            return

        self.mout = self.runfort.share_input_data(self.output_folder)
#        print('self.output_folder',self.output_folder)
        self.show_output()

    def show_output(self):
        out = Output(self.root)
        out.load_structure_file(self.output_folder)


if __name__ == "__main__":
    root = tk.Tk()
    root.title("Chemical Insight & Symmetric VB (CISVB) Structures Generation Software")
    root.geometry("600x920")
    root.configure(background="lightblue")
    GlobVar.configure_styles()
    inputc = Ctrl_Input(root)
    manager = class_manager(root, inputc)

    style_colour_frame = ttk.Style()
    style_colour_frame.configure("Colour_Frame.TFrame", background="lightblue")
    frame1 = ttk.Frame(root, style="Colour_Frame.TFrame", padding="10")
    frame1.grid(row=5, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
    frame2 = ttk.Frame(root, style="Colour_Frame.TFrame", padding="10")
    frame2.grid(row=6, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

    Run_label = ttk.Label(frame1, text='Run Calculation:', style='Colour_Label3.TLabel')
    Run_label.grid(row=1, column=0, padx=15, pady=10)
    Run_button = ttk.Button(frame1, text="RUN", command=lambda: manager.Run_fort_subs())
    Run_button.grid(row=2, column=3, columnspan=3, padx=50, pady=10)
    Run_button.config(state=tk.DISABLED)  # Initially disable the button

    keywd_button = ttk.Button(frame1, text="KEYWDS", command=lambda: manager.create_keywd(Run_button))
    keywd_button.grid(row=0, column=3, padx=15, pady=10)
    keywd_button.config(state=tk.DISABLED)  # Initially disable the button

    orbital_button = ttk.Button(frame1, text="ORBITALS", command=lambda: manager.create_orb(keywd_button))
    orbital_button.grid(row=0, column=4, padx=15, pady=10)
    orbital_button.config(state=tk.DISABLED)  # Initially disable the button
    inputc.orbital_button = orbital_button

    close_label = ttk.Label(frame2, text='Finish Calculation:', style='Colour_Label3.TLabel')
    close_label.grid(row=0, column=1, padx=15, pady=10)
    close_button = ttk.Button(frame2, text="FINISH", command=manager.finish)
    close_button.grid(row=1, column=2, columnspan=3, padx=20, pady=10)

    root.mainloop()
