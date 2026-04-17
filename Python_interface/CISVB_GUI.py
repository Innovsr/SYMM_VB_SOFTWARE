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
import tkinter.font as tkFont
import math
from math import acos, degrees
from itertools import combinations
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D, proj3d
import symm_str
import webbrowser
from collections import defaultdict


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
    all_atoms = []      # List to store atom data
    set_id = 0               # output set id to store present (on show) set number
    symm_key = 0
    total_atoms = None
    CovIon = None
    num_sets = None
    Rum_Ch = None
    IAB_flag = True
    geo_unit = None
    #Gl_atoset = np.zeros((total_atoms, num_electron), dtype=int)
    Gl_atoset = None 
    #Gl_activeatoms = np.zeros(30, dtype=np.int32)
    Gl_activeatoms = None
    active_orbitals = []
    quality_fac = []
    ciflg = ''
    at_list_bold = [            # List of Atoms according to periodic table
        'H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE', 'NA', 'MG', 'AL',
        'SL', 'P', 'S', 'CL', 'AR', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NR',
        'CU', 'ZN', 'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC',
        'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', 'CE',
        'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA',
        'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA'
    ]

    at_covrad = {
            'H': 37, 'HE': 32, 'LI': 134, 'BE': 90, 'B': 82, 'C': 77, 'N': 75, 'O': 73, 
            'F': 71, 'NE': 69,'NA': 154, 'MG': 130, 'AL': 118, 'SL': 111, 'P': 106, 'S': 102,
            'CL': 99, 'AR': 97, 'K': 196, 'CA': 174,'SC': 144, 'TI': 136, 'V': 125, 'CR': 127,
            'MN': 139, 'FE': 125, 'CO': 126, 'NR': 121, 'CU': 138, 'ZN': 131,'GA': 126, 'GE': 122,
            'AS': 119, 'SE': 116, 'BR': 114, 'KR': 110, 'RB': 211, 'SR': 192, 'Y': 162, 'ZR': 148,
            'NB': 137, 'MO': 145, 'TC': 156, 'RU': 126, 'RH': 135, 'PD': 131, 'AG': 153, 'CD': 148,
            'IN': 144, 'SN': 141, 'SB': 138, 'TE': 135, 'I': 133, 'XE': 130, 'CS': 225, 'BA': 198,
            'LA': 169, 'CE': 204, 'PR': 203, 'ND': 201, 'PM': 199, 'SM': 198, 'EU': 198, 'GD': 196,
            'TB': 194, 'DY': 192, 'HO': 192, 'ER': 189, 'TM': 190, 'YB': 187, 'LU': 160, 'HF': 150,
            'TA': 138, 'W': 146, 'RE': 159, 'OS': 128, 'IR': 137, 'PT': 128, 'AU': 144, 'HG': 149,
            'TL': 148, 'PB': 147, 'BI': 146, 'PO': 140, 'AT': 150, 'RN': 145, 'FR': 260, 'RA': 221
            }

    colors = {
            'H': 'white', 'HE': 'cyan', 'LI': 'purple', 'BE': 'darkgreen', 'B': 'salmon',
            'C': 'gray', 'N': 'blue', 'O': 'red', 'F': 'green', 'NE': 'cyan',
            'NA': 'blue', 'MG': 'green', 'AL': 'gray', 'SL': 'orange', 'P': 'orange',
            'S': 'yellow', 'CL': 'green', 'AR': 'cyan', 'K': 'purple', 'CA': 'darkgreen',
            'SC': 'gray', 'TI': 'gray', 'V': 'gray', 'CR': 'gray', 'MN': 'gray',
            'FE': 'orange', 'CO': 'white', 'NR': 'gray', 'CU': 'brown', 'ZN': 'gray',
            'GA': 'gray', 'GE': 'gray', 'AS': 'gray', 'SE': 'orange', 'BR': 'darkred',
            'KR': 'cyan', 'RB': 'purple', 'SR': 'darkgreen', 'Y': 'gray', 'ZR': 'gray',
            'NB': 'gray', 'MO': 'gray', 'TC': 'gray', 'RU': 'gray', 'RH': 'gray',
            'PD': 'gray', 'AG': 'silver', 'CD': 'gray', 'IN': 'gray', 'SN': 'gray',
            'SB': 'gray', 'TE': 'gray', 'I': 'purple', 'XE': 'cyan', 'CS': 'purple',
            'BA': 'darkgreen', 'LA': 'gray', 'CE': 'gray', 'PR': 'gray', 'ND': 'gray',
            'PM': 'gray', 'SM': 'gray', 'EU': 'gray', 'GD': 'gray', 'TB': 'gray',
            'DY': 'gray', 'HO': 'gray', 'ER': 'gray', 'TM': 'gray', 'YB': 'gray',
            'LU': 'gray', 'HF': 'gray', 'TA': 'gray', 'W': 'gray', 'RE': 'gray',
            'OS': 'gray', 'IR': 'gray', 'PT': 'gray', 'AU': 'gold', 'HG': 'silver',
            'TL': 'gray', 'PB': 'gray', 'BI': 'gray', 'PO': 'gray', 'AT': 'gray',
            'RN': 'cyan', 'FR': 'purple', 'RA': 'darkgreen'
            }

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
                "Colour_Label3.TLabel", foreground="darkblue", background="lightblue", font=("Arial", 16, "bold")
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
        self.buttons = {}

        
        group_label_info = {
                "Top_label":("frame","Basic Informations About The Molecular System",1,18),
                "bottom_label":("frame4","Orbital Info & Keywords:",5,16),
                "geometry_label":("frame2","Geometry of the Molecular System",0,12),
                "Active_space_label":("frame3","Define Active Space of the System",0,12)
                }

        for label_name, (frame, text, row, font) in group_label_info.items():
            if frame in frame_info:
                label = ttk.Label(self.frames[frame], text=text, style="Colour_Label3.TLabel",
                                  font=("Airtel", font, "bold"))
                label.grid(row=row, column=0, columnspan=3, padx=25, pady=25)

        # create Labels
        label_info ={
                "brows_geo_label":("frame2","Brows to Upload Geometry", 2, 0, 
                                   "If you have the geometry saved in any dat file you can brows\n"
                                   "that file and insert the geometry here using 'Brows' button; in\n" 
                                   "the geometry file you should have 4 or 5 columns first column: \n"
                                   "Atoms, second column: atomic numbers,third column: x coordinates,\n"
                                   "fourth column: y coordinates, fifth column: z coordinates"),
                "manual_geo_label":("frame2", "Insert Geometry Manually", 3, 0,
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
                "help_button":("frame", "Help", self.Call_help_info, 0, 2,
                               "Please click to get help in using this software"),
                "brows_geo_button":("frame2", "Brows", self.read_geometry, 2, 1,
                                    "If you have the geometry saved in any dat file you can brows\n"
                                    "that file and insert the geometry here using 'Brows' button; in the \n"
                                    "geometry file you should have 4 or 5 columns first column: Atoms, \n"
                                    "second column: atomic numbers, third column: x coordinates, fourth\n" 
                                    "column: y coordinates, fifth column: z coordinates"),
                "manual_geo_button":("frame2", "Geometry", self.insert_geo_manually, 3, 1, 
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
                sticky_opt = 'e' if text == 'Help' else 'w'
                button.grid(row=row, column=column, sticky=sticky_opt, padx=5, pady=5)
                self.balloon.bind(button, tooltip)
                self.buttons[button_name] = button
                if button_name in ('brows_geo_button','manual_geo_button'): 
                    button.config(state=tk.DISABLED)  # Initially disable the button

        self.create_ctrl_pans()

    def Call_help_info(self):
        a=Help_info()
        a.open_help_window()

    def create_geo_unit(self):
        '''
        Its creats the radio button to get the information about the units of geometry provided.
        there are only two options: Bohr or Angstrom.
        '''
        units = ["Bohr", "Angs"]
        label = ttk.Label(self.frames["frame2"], text="Unit of the Geometry Data", style="Colour_Label.TLabel")
        label.grid(row=1, column=0, sticky=tk.W, padx=15, pady=5)
        for i, unit in enumerate(units, start=1):
            button = ttk.Radiobutton(
                self.frames["frame2"],
                text=unit,
                value=unit,
                variable=self.unit_type,
                command=self.update_geo_unit,
                style="Custom.TRadiobutton"
            )
            button.grid(row=1, column=i, padx=10, pady=10)

    def update_geo_unit(self):
        GlobVar.geo_unit = self.unit_type.get()
        if GlobVar.geo_unit:
            self.unit_type_entry = True
            for name in ('brows_geo_button', 'manual_geo_button'):
                if name in self.buttons:
                    self.buttons[name].config(state=tk.NORMAL)
            print('geo_unit**', GlobVar.geo_unit)
            return (GlobVar.geo_unit)

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
            #pattern = re.compile(r"([A-Z][a-z]?)(\d*)\((\d*)([+-])\)")
            #pattern = re.compile(r"([A-Z][a-z]?)(\d*)\((\d*)([+-])\)")
            pattern = re.compile(r"([A-Z][a-z]?)(\d*)(?:\((\d*)([+-])\))?")
            print("pattern",pattern)
            # Remove underscores for handling O_H2 as OH2
            input_string = mol_entry.replace("_", "")
            # Find all matches in the string
            matches = pattern.findall(input_string)
            print('matches',matches)
            charge_num, charge_sign = None, None
            for element, count, num, sign in matches:
                print('element, count, num, sign',element, count, num, sign )
                if element:
                    count = int(count) if count else 1 # If no count is given, assume it's 1
#                    molecule[element] = molecule.get(element, 0) + count
                else:
                    raise ValueError(f"Element {element} is not spelled Correctly.")

                prev_count, prev_charge_num, prev_charge_sign = molecule.get(
                    element, (0, None, None)
                )

                if sign:
                    charge_num = int(num) if num else 1
                    charge_sign = sign
                else:
                    charge_num = prev_charge_num
                    charge_sign = prev_charge_sign

                molecule[element] = (
                    prev_count + count,
                    charge_num,
                    charge_sign
                )

               # molecule[element] = (molecule.get(element, 0) + count, charge_num, charge_sign)

                print("molecule",molecule)

            total = 0
            for element, (count, charge_num, charge_sign) in molecule.items():
                if element in GlobVar.at_list_bold:
                    atomic_number = GlobVar.at_list_bold.index(element) + 1  # Atomic number = index + 1
                    total += atomic_number * count
                    print("total",total)
                if charge_sign == '+':
                    total -= charge_num 
                    print("total+",total)
                if charge_sign == '-':
                    total += charge_num 
                    print("total-",total)
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
                ).grid(row=2, column=2, columnspan=2)

    def insert_geo_manually(self):
        """Allows manual insertion of geometry data via Read_Geo."""
        GlobVar.readgeo = Read_Geo(self.file_path)  # Initialize the Read_Geo class
        GlobVar.readgeo.insert_geo(self.root)  # Call insert_geo method
        ttk.Button(
            self.frames["frame2"],
            text="View Geometry",
            command=GlobVar.readgeo.display_geometry
            ).grid(row=3, column=2, columnspan=2)

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
        self.symat = []      # list to store atom names
        self.symatno = []    # list to store atomic numbers
        self.coordx = []     # list to store atoms x coordinates
        self.coordy = []     # list to store atoms y coordinates
        self.coordz = []     # list to store atoms z coordinates
        self.atoms = []

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
        bohrtoangs = 1.0
        if GlobVar.geo_unit == 'Bohr':
            bohrtoangs = 0.529177

        try:
            with open(self.file_path, 'r') as file:
                lines = file.readlines()
                nline = 0
            for line in lines:
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                columns = line.split()
                nline += 1
                if len(columns) == 4:
                    ''' 4-column format: atom name, x, y, z'''
                    atom_name = columns[0].capitalize()  # Capitalize atom name
                    x = self.Deconvert(columns[1])*bohrtoangs
                    y = self.Deconvert(columns[2])*bohrtoangs
                    z = self.Deconvert(columns[3])*bohrtoangs
                    self.atoms.append({
                        "sl_num":nline,
                        "atom": atom_name,
                        "x": x,
                        "y": y,
                        "z": z
                    })
                elif len(columns) == 5:
                    ''' 5-column format: atom name, atomic number, x, y, z'''
                    atom_name = columns[0].capitalize()  # Capitalize atom name
                    atomic_number = self.Deconvert(columns[1])
                    x = self.Deconvert(columns[2])*bohrtoangs
                    y = self.Deconvert(columns[3])*bohrtoangs
                    z = self.Deconvert(columns[4])*bohrtoangs
                    self.atoms.append({
                        "sl_num":nline,
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

            print('atom', self.atoms)
            GlobVar.total_atoms = len(self.atoms)
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
            GlobVar.all_atoms = self.atoms
            return self.symat, self.coordx, self.coordy, self.coordz, self.symatno
        except Exception as e:
            print(f"An error occurred: {e}")
            return None
        return (self.atoms)

    def insert_geo(self, root):
        geo_window = tk.Toplevel(root)
        geo_window.title("Geometry")
        geo_window.geometry("650x650")
        geo_window.configure(background="lightblue")
        self.entry_widgets = []
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
        self.frame_entry = ttk.Frame(canvas, style="Colour_Frame.TFrame", padding=10)

        # Configure scrollable area
        self.frame_entry.bind(
            "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=self.frame_entry, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        label = ttk.Label(frame1, text="Please insert the number of atoms: ", style="Colour_Label.TLabel")
        label.grid(row=0, column=0, padx=10, pady=10, sticky=tk.W)

        atom_num_entry = ttk.Entry(frame1, width=10)
        atom_num_entry.grid(row=0, column=1, padx=10, pady=10)

        button_Enter = ttk.Button(frame1, text="Enter", command=lambda:self.create_pan(atom_num_entry))
        button_Enter.grid(row=0, column=3, padx=10, pady=10, sticky=tk.W)
        button_Insert = ttk.Button(frame3, text="Insert", command=self.fetch_data)
        button_Insert.grid(row=0, column=3, padx=10, pady=10, sticky=tk.W)
        button_Close = ttk.Button(frame3, text="Close", command=geo_window.destroy)
        button_Close.grid(row=1, column=3, padx=10, pady=10, sticky=tk.W)
        return (self.atoms)

    def create_pan(self, atom_num_entry):
        self.atom_num = int(atom_num_entry.get())

        for widget in self.frame_entry.winfo_children():
            widget.destroy()
        self.entry_widgets.clear()

        label1 = ttk.Label(self.frame_entry, text=" Atom ", background="lightblue")
        label1.grid(row=1, column=0, padx=10)
        label2 = ttk.Label(self.frame_entry, text=" X Coordinates ", background="lightblue")
        label2.grid(row=1, column=1, padx=10)
        label3 = ttk.Label(self.frame_entry, text=" Y Coordinates ", background="lightblue")
        label3.grid(row=1, column=2, padx=10)
        label4 = ttk.Label(self.frame_entry, text=" Z Coordinates ", background="lightblue")
        label4.grid(row=1, column=3, padx=10)

        for i in range(self.atom_num):
            Atom = ttk.Entry(self.frame_entry, width=10)
            Atom.grid(row=2+i, column=0, padx=10, pady=10)

            X_coord = ttk.Entry(self.frame_entry, width=15)
            X_coord.grid(row=2+i, column=1, padx=10, pady=10)

            Y_coord = ttk.Entry(self.frame_entry, width=15)
            Y_coord.grid(row=2+i, column=2, padx=10, pady=10)

            Z_coord = ttk.Entry(self.frame_entry, width=15)
            Z_coord.grid(row=2+i, column=3, padx=10, pady=10)

            self.entry_widgets.append((Atom, X_coord, Y_coord, Z_coord))

    def fetch_data(self):
        self.atoms.clear()
        self.symat.clear()
        self.coordx.clear()
        self.coordy.clear()
        self.coordz.clear()
        self.symatno.clear()

        bohrtoangs = 1.0
        if GlobVar.geo_unit == 'Bohr':
            bohrtoangs = 0.529177

        for atom, x_entry, y_entry, z_entry in self.entry_widgets:
            atom_name = atom.get().strip()
            if atom:  # Only process if the atom field is not empty
                try:
                    x = float(self.Deconvert(x_entry.get().strip()))*bohrtoangs
                    y = float(self.Deconvert(y_entry.get().strip()))*bohrtoangs
                    z = float(self.Deconvert(z_entry.get().strip()))*bohrtoangs
                except ValueError:
                    messagebox.showerror("Input Error", "Coordinates must be valid numbers.")
                    return

                self.atoms.append({
                    "atom": atom_name,
                    "x": x,
                    "y": y,
                    "z": z
                })

        GlobVar.total_atoms = self.atom_num
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
                messagebox.showerror("Name Error",f"Element {element} is not spelled Correctly.\n"
                                    "or, maybe its not in Uppercase or Uppercase-Lowercase letters "
                                    "combination or maybe it's above 88 elements of the periodic table \n"
                                    "we only consider fistr 88 elements of the periodic table")
                return

        #print("self.symat, self.coordx, self.coordy, self.coordz, self.symatno", self.symat, self.coordx,\
        #self.coordy, self.coordz, self.symatno)
        GlobVar.geometry_inserted = True
        GlobVar.all_atoms = self.atoms

        return self.symat, self.coordx, self.coordy, self.coordz, self.symatno


    def display_geometry(self):
        ''' Create a Toplevel window to display the file geometry '''
        if GlobVar.geometry_inserted is True:
            geo_file_display = tk.Toplevel()
            if (GlobVar.molecule_string):
                geo_file_display.title(f"{GlobVar.molecule_string}-geometry")
            else:
                geo_file_display.title("-Geometry-")
            geo_file_display.geometry("500x300")  # Optional: Set the size of the Toplevel window

            geo_file_display_frame = tk.Frame(geo_file_display)
            geo_file_display_frame.grid(row=0,column=0, sticky='nsew')
            geo_file_display.grid_rowconfigure(0, weight=1)
            geo_file_display.grid_columnconfigure(0, weight=1)
            self.geo_mol_show_frame = tk.Frame(geo_file_display)
            self.geo_mol_show_frame.grid(row=1,column=0, sticky='nsew')
            close_frame = tk.Frame(geo_file_display)
            close_frame.grid(padx=10, pady=10, sticky='nsew')

            style = ttk.Style()
            style.configure("Treeview", font=("Helvetica", 12))
            # Define the Treeview widget with columns
            columns = ("atom", "x", "y", "z")
            tree = ttk.Treeview(geo_file_display, columns=columns, show="headings")
            #tree.pack(expand=True, fill="both")
            tree.grid(row=0, column=0, sticky="nsew")

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

        # Add a scrollbar for better usability
        scrollbar = ttk.Scrollbar(geo_file_display, orient="vertical", command=tree.yview)
        tree.configure(yscrollcommand=scrollbar.set)
        tree.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")


        geo_file_display.grid_rowconfigure(0, weight=1)
        geo_file_display.grid_columnconfigure(0, weight=1)
        self.display_molecule()

        close_button = tk.Button(close_frame, text='Close', command=geo_file_display.destroy)
        close_button.pack(expand=True)

    def display_molecule(self):
        if GlobVar.geometry_inserted is True:
            atoms_new = []
            atoms_new = self.atoms
            mol_display = tk.Toplevel()
            if (GlobVar.molecule_string):
                mol_display.title(f"{GlobVar.molecule_string}")
            else:
                mol_display.title("-Molecule-")
            mol_display.geometry("700x700")
            mol_display.configure(background="lightblue")

            #self.atoms = GlobVar.atoms

            mol_frame = tk.Frame(mol_display)
            mol_frame.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

            button_frame = tk.Frame(mol_frame)
            button_frame.pack(fill=tk.BOTH, expand=True)

            center_frame = tk.Frame(button_frame)
            center_frame.pack(expand=True)

            canvas_frame = tk.Frame(mol_frame)
            canvas_frame.pack(fill=tk.BOTH, expand=True)

            box_frame = tk.Frame(mol_frame, width=100, height=100, bg='blue')
            box_frame.pack(pady=10)
            box_frame.pack_propagate(False)
            
            fig = plt.Figure(figsize=(6, 6))
            fig.patch.set_facecolor('lightblue')
            ax = fig.add_subplot(111, projection='3d')
            ax.set_facecolor('lightblue')
            
            # drawing axis
            ax.plot([0, 2], [0, 0], [0, 0], color='blue', linewidth=2, linestyle=':')
            ax.plot([0, 0], [0, 2], [0, 0], color='blue', linewidth=2, linestyle=':')
            ax.plot([0, 0], [0, 0], [0, 2], color='blue', linewidth=2, linestyle=':')
            ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio

            # drawing atoms
            old_colors = []
            self.selected_atoms = []
            atom_scatters = []
            scale_factor = 0.5
            for atom1 in atoms_new:
                x1 = atom1['x']
                y1 = atom1['y']
                z1 = atom1['z']
                symbol = atom1['atom']

                radius = GlobVar.at_covrad.get(symbol)  # fallback value if element missing
                size = (radius * scale_factor)**2
                print('size of the atom',size)

                scatters = ax.scatter(x1, y1, z1,
                           color=GlobVar.colors.get(atom1['atom'], 'green'),
                           s=size)

                atom_scatters.append(scatters)
                ax.text(x1, y1, z1 + 0.3, atom1['atom'], color='black', fontsize=10, ha='center', va='bottom')
                clr = GlobVar.colors.get(atom1['atom'])
                old_colors.append(clr)
            print('scatters',scatters)

            for atom1, atom2 in combinations(atoms_new, 2):
                x1, y1, z1 = atom1['x'], atom1['y'], atom1['z']
                x2, y2, z2 = atom2['x'], atom2['y'], atom2['z']

                distance = math.sqrt(
                        (x2-x1)**2 +
                        (y2-y1)**2 +
                        (z2-z1)**2
                        )

                r1 = GlobVar.at_covrad.get(atom1['atom'])
                r2 = GlobVar.at_covrad.get(atom2['atom'])

                if r1 is not None and r2 is not None:
                    covrad = (r1 + r2)/100.0
                    if distance <= covrad:
                        ax.plot([x1, x2], [y1, y2], [z1, z2], color='black', linewidth=2.5)
                else:
                    messagebox.showerror("Not Found",f"covalent radious of {atom1} and/or {atom2}\n"
                                         "are not found in the list")
                            
            
            # Make panes transparent (no background)
            ax.xaxis.set_pane_color((1, 1, 1, 0))
            ax.yaxis.set_pane_color((1, 1, 1, 0))
            ax.zaxis.set_pane_color((1, 1, 1, 0))
            
            # Hide grid lines
            ax.grid(False)
            
            # Hide tick marks and labels
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            
            # Hide box edges (axis lines)
            try:
                ax.w_xaxis.line.set_color((0, 0, 0, 0))
                ax.w_yaxis.line.set_color((0, 0, 0, 0))
                ax.w_zaxis.line.set_color((0, 0, 0, 0))
            except AttributeError:
                # For newer Matplotlib versions
                for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
                    axis.line.set_alpha(0)

            # Add axis labels
            ax.text(2, 0, 0, 'X', color='blue', fontsize=14, fontweight='bold')
            ax.text(0, 2, 0, 'Y', color='blue', fontsize=14, fontweight='bold')
            ax.text(0, 0, 2, 'Z', color='blue', fontsize=14, fontweight='bold')

            
            canvas = FigureCanvasTkAgg(fig, master=canvas_frame)
            canvas_widget = canvas.get_tk_widget()
            canvas_widget.pack(fill=tk.BOTH, expand=True)

            def on_click(event):
                if event.inaxes != ax:
                    return

                # Get 2D click position
                x2, y2 = event.x, event.y

                # Get projection of 3D points
                coords_2d = []
                for atom in atoms_new:
                    x3, y3, z3 = atom['x'], atom['y'], atom['z']
                    x_proj, y_proj, _ = proj3d.proj_transform(x3, y3, z3, ax.get_proj())
                    x_disp, y_disp = ax.transData.transform((x_proj, y_proj))
                    coords_2d.append((x_disp, y_disp))

                coords_2d = np.array(coords_2d)
                distances = np.linalg.norm(coords_2d - np.array([x2, y2]), axis=1)
                dist = min(distances)
                index = np.argmin(distances)
                atom_clicked = atoms_new[index]
                r1 = GlobVar.at_covrad.get(atom_clicked['atom'])
                print('dist, index',dist, index)
                if dist <= math.sqrt(size):
                    atom_scatters[index].set_color('magenta')
                    fig.canvas.draw_idle()
                    self.selected_atoms.append(atom_clicked)
                    #print('selected_atoms',selected_atoms)

                return (self.selected_atoms)

                # Connect the event
            fig.canvas.mpl_connect('button_press_event', on_click)
            print('selected_atoms',self.selected_atoms)
            
            canvas.draw()

            length = tk.Button(center_frame, text='length',
                               command=lambda:self.calculate_length(
                                                               atom_scatters,
                                                               old_colors, fig)
                               )
            length.pack(side='left', padx=10)
            ang = tk.Button(center_frame, text='ang',
                               command=lambda:self.calculate_angle(
                                                               atom_scatters,
                                                               old_colors, fig)
                            )                        
            ang.pack(side='left', padx=10)
            dang = tk.Button(center_frame, text='di-ang',
                               command=lambda:self.calculate_diangle(
                                                               atom_scatters,
                                                               old_colors, fig)
                            )                        
            dang.pack(side='left', padx=10)

            close_btn = tk.Button(box_frame, text="Close", command=mol_display.destroy)
            close_btn.pack(expand=True)

    def calculate_length(self, atom_scatters, old_colors, fig):
        if len(self.selected_atoms) != 2:
            messagebox.showerror('Selection Error','Please select any two atoms to calculate Bond Length')
            return

        for atom1, atom2 in combinations(self.selected_atoms, 2):
            dx = atom2['x'] - atom1['x']
            dy = atom2['y'] - atom1['y']
            dz = atom2['z'] - atom1['z']

            distance = math.sqrt(dx**2 + dy**2 + dz**2)

        messagebox.showinfo("Bond Length:",f"Distance between {atom1['atom']} and\n"
                            f"{atom2['atom']}: {distance:.4f} Å")

        for i, scatter in enumerate(atom_scatters):
            scatter.set_color(old_colors[i])  # Restore original
        fig.canvas.draw_idle()
        self.selected_atoms = []


    def calculate_angle(self, atom_scatters, old_colors, fig):
        if len(self.selected_atoms) != 3:
            messagebox.showerror('Selection Error','Please select any three atoms to calculate Bond Angle')
            return
        for atom1, atom2, atom3 in combinations(self.selected_atoms, 3):
            a = np.array([atom1['x'], atom1['y'], atom1['z']])
            b = np.array([atom2['x'], atom2['y'], atom2['z']])
            c = np.array([atom3['x'], atom3['y'], atom3['z']])

        # the vectpr formula is used here : cos(θ)= u.v/|u||v|
        ba = np.array(a) - np.array(b) # ba = u
        bc = np.array(c) - np.array(b) # bc = v
        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle_rad = acos(np.clip(cosine_angle, -1.0, 1.0))  # Avoid rounding errors
        angle_deg = degrees(angle_rad)
        messagebox.showinfo("Mesured Angle:",f"Angle within {atom1['atom']},\n"
                            f"{atom2['atom']} and {atom3['atom']}: {angle_deg} D")

        for i, scatter in enumerate(atom_scatters):
            scatter.set_color(old_colors[i])  # Restore original
        fig.canvas.draw_idle()
        
        self.selected_atoms = []
        #return selected_atoms

    def calculate_diangle(self, atom_scatters, old_colors, fig):
        if len(self.selected_atoms) != 4:
            messagebox.showerror('Selection Error','Please select any four atoms to calculate Dihedral Angle')
            return
        for atom1, atom2, atom3, atom4 in combinations(self.selected_atoms, 4):
            a = np.array([atom1['x'], atom1['y'], atom1['z']])
            b = np.array([atom2['x'], atom2['y'], atom2['z']])
            c = np.array([atom3['x'], atom3['y'], atom3['z']])
            d = np.array([atom3['x'], atom3['y'], atom3['z']])

        b0 = -1.0 * (np.array(b) - np.array(a))
        b1 = np.array(c) - np.array(b)
        b2 = np.array(d) - np.array(c)

        b1 /= np.linalg.norm(b1)

        v = b0 - np.dot(b0, b1) * b1
        w = b2 - np.dot(b2, b1) * b1

        x = np.dot(v, w)
        y = np.dot(np.cross(b1, v), w)
        dhang = degrees(np.arctan2(y, x))
        messagebox.showinfo("Mesured Dihedral Angle:",f"Angle within {atom1['atom']},\n"
                            f"{atom2['atom']}, {atom3['atom']} and {atom4['atom']}: {dhang} D")

        for i, scatter in enumerate(atom_scatters):
            scatter.set_color(old_colors[i])  # Restore original
        fig.canvas.draw_idle()
        
        self.selected_atoms = []

    def get_geometry_data(self):
        return self.symat, self.coordx, self.coordy, self.coordz, self.symatno


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
        self.atoset = np.zeros((GlobVar.total_atoms, GlobVar.num_electron), dtype=np.int32)
        self.norbsym_py = np.empty(4, dtype=np.int32)
        self.orbsym = np.zeros((GlobVar.num_electron, GlobVar.num_electron), dtype=np.int32)
        self.activeatoms = np.zeros(GlobVar.total_atoms, dtype=np.int32)
        self.atn_vector = np.zeros(GlobVar.total_atoms, dtype=np.int32)
        self.description_orb = ""
        Pmw.initialise(root)
        self.balloon = Pmw.Balloon(root)

        self.create_orbital_section()
        #Gl_activeatoms = np.zeros(30, dtype=np.int32)
        GlobVar.Gl_activeatoms = self.activeatoms

    def create_orbital_section(self):
        if GlobVar.orbital_input is False:
            self.orbital_frame = tk.Toplevel(root, padx=5, pady=5)
            self.orbital_frame.title(f"{GlobVar.molecule_string} orbital info")
            self.orbital_frame.geometry("560x560")
            self.orbital_frame.configure(background="lightblue")
            GlobVar.orbital_input = True
            atoset_created = False

            for widget in self.orbital_frame.winfo_children():
                widget.destroy()

            top_label = ttk.Label(
                    self.orbital_frame, text=f"Number of active orbitals: {self.num_orbital}",
                    style="Colour_Label.TLabel"
                    )
            top_label.grid(row=0, column=0, columnspan=3, pady=15)

            label_info={
                    'label1':("Active Orbital",0,"Active orbitals are listed below"),
                    'label2':("Atom Number",1,"Put the active atom number as given in the geometry file\n"
                              "associated with the each active orbitals"),
                    'label3':("Orbital Type",2,"Put the active orbital types; It can be given\n"
                              "as: s or sig, px... or pi1, py... or pi2, pz... or pi3")
                }

            for label, (text, col, tooltip) in label_info.items():
                label = ttk.Label(self.orbital_frame, text=text, style='Colour_Label1.TLabel')
                label.grid(row=1, column=col, padx=5, sticky='ew')
                self.balloon.bind(label, tooltip)

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
                    self.orbital_frame, text="Insert", command=self.validate_and_store_orbital_data)
            insert_button.grid(row=3, column=0, columnspan=2, pady=10)

            self.show_button = ttk.Button(
                    self.orbital_frame, text="View Orbs", command=self.call_drawing_molecule)
            self.show_button.grid(row=3, column=1, columnspan=2, pady=10)

            if not atoset_created:
                self.show_button.config(state=tk.DISABLED)  # Initially disable the button

            close_button = ttk.Button(self.orbital_frame, text="Close", command=self.destroy_orbs)
            close_button.grid(row=4, column=0, pady=10, columnspan=3)

            for i in range(self.num_orbital):
                label4 = ttk.Label(
                        frame1, text=f"active orbital {i+1} ", style="Colour_Label.TLabel")
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
#            self.assoatm_entry.insert(0, atom_number)
#            self.atm_entry.append(self.assoatm_entry)

    def call_drawing_molecule(self):
        d=drawing_molecule()
        d.VB_view()

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
            if i < len(GlobVar.orbital_data):
                GlobVar.orbital_data[i]["atom_number"] = atom_number
                GlobVar.orbital_data[i]["orbital_type"] = orbital_type
            else:
                GlobVar.orbital_data.append({
                    "atom_number": atom_number,
                    "orbital_type": orbital_type
                })

        #messagebox.showinfo("Success", "Orbital inputs validated and stored successfully.")
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
            GlobVar.active_orbitals.append(i+GlobVar.num_iao)
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
            elif 'sig' in data["orbital_type"].lower():
                sig_type = sig_type + 1
                sorbs.append(i+GlobVar.num_iao)
            else:
                messagebox.showerror("Input Error", f"The type of the orbital {i} is unknown. \n"
                                     "Please put s, px, py or pz otherwise put sig or sigma,\n"
                                     " pi1, pi2, pi3 if the direction of pi orbs are different"
                                     )

        if all(x != 0 for x in [py_type, px_type, pz_type]):
            messagebox.showerror("Input Error", "You have put three different pi directions \n"
                                 "among them one should be Sigma, Please check and rectify"
                                 )


        #norbsym_py = np.zeros(4, dtype=np.int32)
        #norbsym = [
        #        px_type,
        #        py_type,
        #        pz_type,
        #        sig_type
        #        ]

        #self.norbsym_py = np.zeros(10, dtype=int)
        #self.norbsym_py[:len(norbsym)] = norbsym
        #self.norbsym_py = norbsym
        self.norbsym_py = np.array(
            [px_type, py_type, pz_type, sig_type],
            dtype=np.int32,
            order='F'
        )

        max_len = GlobVar.num_electron

        # Pad lists to be of equal length
        pxorbs_padded = pxorbs + [0] * (max_len - len(pxorbs))
        pyorbs_padded = pyorbs + [0] * (max_len - len(pyorbs))
        pzorbs_padded = pzorbs + [0] * (max_len - len(pzorbs))
        sorbs_padded = sorbs + [0] * (max_len - len(sorbs))
        self.orbsym[0, :] = pxorbs_padded
        self.orbsym[1, :] = pyorbs_padded
        self.orbsym[2, :] = pzorbs_padded
        self.orbsym[3, :] = sorbs_padded

        print('norbsym', self.norbsym_py)
#        print('sig_type, px_type, py_type, pz_type', sig_type, px_type, py_type, pz_type )

        # Checking how many types are non-zero
        GlobVar.type_orb_count = sum(1 for count in [sig_type, px_type, py_type, pz_type] if count > 0)

        self.atoset = self.create_matrix()
        atoset_created = True
        """atoset is matrix where each row represents an atom according to the geometry, if the atom ia
        an active atom, the corresponding column get '1' otherwise '0'. and the next columns contain the
        corresponding active orbital numbers associated with that atom."""
#        for i 
#        GlobVar.active_atom_coords
        self.keywd_button.config(state=tk.NORMAL)
        self.show_button.config(state = tk.NORMAL)
        GlobVar.Gl_atoset = self.atoset
        print('atoset',self.atoset)
#        self.orbital_button.config(state=tk.DISABLED)  # Initially disable the button

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
        matrix = np.zeros((GlobVar.total_atoms, GlobVar.num_electron), dtype=int)

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
        for attn in self.atn_vector:
            if attn > 1:
                GlobVar.IAB_flag = False
                break

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
        self.get_set_order_window = None
        self.frame = None
        self.method_type = tk.StringVar(value="Chemical Insight Set")
        self.ChemInst_set_type = tk.StringVar(value="Single Set")
        self.rumer_set_type = tk.StringVar(value="Single Rumer Set")
        self.str_type = tk.StringVar(value="Covalent")
        self.cheminst_str_type = tk.StringVar(value="Symmetric")
        self.symmetry_set_order_type = tk.StringVar(value="Quality-Arrange")
        self.ovlp_cal = tk.StringVar(value="No")
        self.i=0

        if GlobVar.IAB_flag is True:
            self.IAB_type = tk.StringVar(value="None")
        else:
            self.IAB_type = tk.StringVar(value="1")

        if GlobVar.total_atoms < 3:
            self.NAB_type = tk.StringVar(value="None")
        else:
            self.NAB_type = tk.StringVar(value="2")

        if self.type_orb_count == 1:
            self.SBB_type = tk.StringVar(value="None")
        else:
            self.SBB_type = tk.StringVar(value="3")

        self.PDB_type = tk.StringVar(value="None")
        self.PDR_type = tk.StringVar(value="None")
        self.ChemInst_set_type_entry = True
        #self.cheminst_str_type_entry = False
        self.symmetric_set_order_type_entry = False
        self.prio_structure_entries = []
        self.prio_bond_entries = []
        self.prio_rads_entries = []
        self.prio_lnp_entries = []
        self.PDR_buttons = []
        self.SBB_buttons = []
        self.PDB_data = []
        Pmw.initialise(self.root)
        self.balloon = Pmw.Balloon(self.root)
        self.mout_number = 1
        self.bond_number = 0
        self.pref_rad_lp_py = []
        self.pref_lnp_py = []
        self.numlp_py = []
        self.numradlp_py = []
        self.pref_rad_py = []
        self.numudr_py = 0
        self.numrad_py = 0 
        self.radio_widgets = []

    def create_keywd_pane(self):
        if self.keywd_window is None:
            self.keywd_window = tk.Toplevel(self.root, padx=10, pady=10)
            self.keywd_window.title("Spatial Keyword Inputs")
            self.keywd_window.geometry("860x950")
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
                "VB-set_type_frame": 1,
                #"rumer_set_type_frame": 1,
                "str_type_frame": 2,
                "ChemInst_set_type_frame": 3,
                #"cheminst_str_type_frame": 4,
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
            col = 0
            for frame_name, row in frames_info.items():
                if frame_name == 'end_frame':
                    col = 2
                frame = ttk.Frame(self.keywd_window, style="Colour_Frame.TFrame")
                frame.grid(row=row, column=col, columnspan=5, sticky=tk.W)
                self.frames[frame_name] = frame
           
            labels_info = {
                "method_type_label": ("method_type_frame", "Select One from the Four Main Operations below: ", 0, 6,
                                      "Select any one from the four methods provided.\n"
                                      "Select 'Rumer' to generate a Rumer structure set,\n"
                                      "no need to provide any spatial keywords. For all\n"
                                      "other types of sets please select Chem Inst."),
    
                "empty_label1": ("method_type_frame", "", 2, 4,
                                         ""),

                "keywd_label": ("VB-set_type_frame", " Select Apropreate Keywords Below:", 0, 4,
                                         ""),

                "str_type_label": ("str_type_frame", "Select Type of Valence Bonds: ", 0, 1,
                                   "Select any one from the three types. Select 'Covalent'\n"
                                   "to generate covalent sets only. Select 'Ionic' to generate \n"
                                   "ionic structure sets only. Select both to generate both \n"
                                   "types of structure sets."),

                "ChemInst_set_type_label": ("ChemInst_set_type_frame", "Number of Complete VB sructure Set:", 0, 1,
                                            "Select any one from the four types of Chemical\n"
                                            "insight output sets. Select 'Single Set' if you want the best\n"
                                            "chemically insightful one set, select 'All Best Sets' if you\n"
                                            "want all possible same best quality sets if available, select\n"
                                            "'All Sets' if you want all possible sets but it can be huge\n"
                                            "(adjust maximum output file accordingly). Select 'Eq Bond' if\n"
                                            "you want only equally distributed bond sets, i.e., lowest \n"
                                            "overlapped sets but less chemically meaningful. These sets \n"
                                            "could help read off negative Coulson-Chirgwin weights."),

                "ChemInst_str_type_label": ("VB-set_type_frame", "VB-Set Type:", 1, 1,
                                            "Select any one from the three types of calculations. Select\n"
                                            "'Symmetry' if you want to have symmetric sets. Select \n"
                                            "'Asymetric' if you want to have chemical insight sets. \n"
                                            "Select 'Checksymm' if you want to check some sets are \n"
                                            "symmetric or not"),


                "ovlp_cal_label": ("ovlp_cal_frame", "Estimate Overlap", 0, 1,
                                   "To estimate overlap among the structures \n"
                                   "in each set click yes. The default is set to 'No'"),

                "empty_label2": ("ovlp_cal_frame", "", 1, 4,
                                         ""),

                "prio_label": ("prio_label_frame", "Select Priorities of the Chemical Qualities below:", 0, 5,
                               "The sequence of the default priority is IAB > NAB > SBB and \n"
                               "PBU & PRU are not taken in the calculation; To change this default \n"
                               "priorities please click the button and selevct the numbers for each "),

                "IAB_label": ("diff_qualities_frame", "IAB  PRIORITY: ", 0, 1,
                              "Intra Atomic Bond Priority. Default is set to '1'\n"
                              "If it shows as disabled, it means that the active atoms \n"
                              "of your system contains only one active orbital."),

                "NAB_label": ("diff_qualities_frame", "NAB  PRIORITY: ", 1, 1,
                              "Near Atomic Bond Priority. Default is set to '2'\n"
                              "If it shows as disabled, it means that your active \n"
                              "space contains only two active atoms."),

                "SBB_label": ("diff_qualities_frame", "SBB  PRIORITY: ", 2, 1,
                              "Symmetry Breaking Bond Priority. Default is set to '3'. \n"
                              "If it shows as disabled, it means that the active space \n"
                              "of your system contains only one type of symmetric orbital."),

                "PDB_label": ("diff_qualities_frame", "PDB  PRIORITY: ", 3, 1,
                              "Pre-defined Bond Priority. Default is set to \n"
                              "'None': Not taken into the calculation"),

                "PDR_label": ("diff_qualities_frame", "PDR  PRIORITY: ", 4, 1,
                              "Pre-defined Radical Priority. Default is set to \n"
                              "'None': Not taken into the calculation"),

                "empty_label3": ("mout_frame", "", 1, 4,
                                         ""),
                }
                
            for label_name, (frame_name, text, row, colspn, tooltip) in labels_info.items():
                styles = "Colour_Label.TLabel"
                if label_name in ("IAB_label", "NAB_label", "SBB_label", "PDB_label", "PDR_label"):
                    styles = "Colour_Label2.TLabel"
                if label_name in (
                        "method_type_label",
                        "prio_label",
                        "empty_label1",
                        "empty_label2", 
                        "empty_label3", 
                        "keywd_label"
                        ):
                    styles = "Colour_Label3.TLabel"

                if frame_name in self.frames:
                    label = ttk.Label(self.frames[frame_name], text=text, style=styles)
                    label.grid(row=row, column=0, columnspan=colspn, sticky=tk.W, padx=10, pady=10)
                    self.balloon.bind(label, tooltip)
                    self.labels[label_name] = label

            self.radio_buttons_info = {
                    #"Rumer_Set_Type":["Single Rumer Set", "All Rumer Sets"],
                    "tip_method":["Chemical Insight Set", "Rumer Set", "Equal Bond-Distributed Set", "Check Set if Symmetric"],
                    "str_type":["Covalent", "Ionic", "Both"],
                    "cheminstsets_type":["Single Set", "All Sets", "Best Sets"],
                    "Cheminst_str_type":["Symmetric", "Asymmetric"],
                    "ovlp_types":["Yes", "No"],
                    "IAB_option":["1", "2", "3", "4", "5", "None"],
                    "NAB_option":["1", "2", "3", "4", "5", "None"],
                    "SBB_option":["1", "2", "3", "4", "5", "None"],
                    "PDB_option":["1", "2", "3", "4", "5", "None"],
                    "PDR_option":["1", "2", "3", "4", "5", "None"]
                    }

            self.radio_buttons_variable = {
                    #"Rumer_Set_Type":["rumer_set_type_frame", self.rumer_set_type, 
                    #                  self.Update_Rumer_Set_Type, 0, 1],
                    "tip_method":["method_type_frame", self.method_type, self.update_method_type, 1, 0],
                    "str_type":["str_type_frame", self.str_type, self.Update_Str_Type, 0, 1],
                    "cheminstsets_type":["ChemInst_set_type_frame", self.ChemInst_set_type, 
                                         self.ChemInst_set_type_read, 0, 1],
                    "Cheminst_str_type":["VB-set_type_frame", self.cheminst_str_type, 
                                         self.cheminst_str_type_read, 1, 1],
                    "ovlp_types":["ovlp_cal_frame", self.ovlp_cal, self.get_str_ovlp_keywd, 0, 1],
                    "IAB_option":["diff_qualities_frame", self.IAB_type, self.get_IAB_Priority, 0, 1],
                    "NAB_option":["diff_qualities_frame", self.NAB_type, self.get_NAB_Priority, 1, 1],
                    "SBB_option":["diff_qualities_frame", self.SBB_type, self.get_SBB_Priority, 2, 1],
                    "PDB_option":["diff_qualities_frame", self.PDB_type, self.get_PDB_Priority, 3, 1],
                    "PDR_option":["diff_qualities_frame", self.PDR_type, self.get_PDR_Priority, 4, 1]
                    }

            for button, (frame, variable, callback, row, colstart) in self.radio_buttons_variable.items():
                if button in self.radio_buttons_info:
                    self.create_radiobuttons(
                            self.frames[frame], self.radio_buttons_info[button], variable, 
                            callback, row, colstart )

            button_info = {
                    "PDB_button":["PDB_PDR_Button_frame", " Insert PDB ", self.Insert_PDB, 0, 0,
                                  "Click the button to insert your preffered bonds"],
                    "PDR_button":["PDB_PDR_Button_frame", " Insert PDR ", self.Insert_PDR, 0, 1,
                                  "Click the button to insert your preffered bonds"],
                    "prio_struc_button":["priority_str_frame", "Insert Preferred Structure", 
                                         self.Insert_Priority_Str, 0, 0, 
                                         "If you want some structures must be present in the set \n"
                                         "you can insert them by clicking this Button;\n"
                                         "The provided structures must be linearly independent \n"
                                         "otherwise they will not be present in the set together"],
                    "insert_button":["end_frame", "Insert", self.get_keywds, 0, 3,
                                     "Click the button to insert all values"],
                    "close_button":["end_frame", "Close", self.keywd_window.destroy, 1, 3,
                                     "Click the button to closr this Keyword pane"],
                    "more_button":["VB-set_type_frame", "More", self.get_set_order, 1, 4,
                                     "Click for more options"],
                    "more_button1":["ChemInst_set_type_frame", "More", self.get_mout_number, 0, 5,
                                     "Click for more options"]
                    }
            for button_name, (frame, text, command, row, column, tooltip) in button_info.items():
                button = ttk.Button(self.frames[frame], text=text, command=command)
                button.grid(row=row, column=column, padx=10, pady=10, sticky=tk.W)
                self.balloon.bind(button, tooltip)
                self.button[button_name] = button

                if button_name in ("PDB_button", "PDR_button"):
                    button.config(state=tk.DISABLED)  # Initially disable the button


    def create_radiobuttons(self, frame, options, variable, callback, row, colstart):
        ''' Create radio buttons; check conditions and disabled specific buttons accordingly '''
        for i, option in enumerate(options, start=colstart):
            button = ttk.Radiobutton(
                frame,
                text=option,
                value=option,
                variable=variable,
                command=callback,
                style="Custom.TRadiobutton"
            )
            button.grid(row=row, column=i, padx=10, pady=10)
            self.radio_widgets[option] = button
            #self.disable_widgets_in_frame(self.frames["rumer_set_type_frame"])
            if GlobVar.multiplicity == 1:
                if options is self.radio_buttons_info["PDR_option"]:
                    button.config(state=tk.DISABLED)
            if self.type_orb_count == 1:
                if options is self.radio_buttons_info["SBB_option"]:
                    button.config(state=tk.DISABLED)
            if GlobVar.total_atoms < 3:
                if options is self.radio_buttons_info["NAB_option"]:
                    button.config(state=tk.DISABLED)
            if GlobVar.IAB_flag is True:
                if options is self.radio_buttons_info["IAB_option"]:
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
        if PDB_Priority is None:
            self.button["PDB_button"].config(state=tk.DISABLED)
        return (PDB_Priority)

    def get_PDR_Priority(self):
        PDR_Priority = self.PDR_type.get()
        print('PDR_Priority', PDR_Priority)
        if PDR_Priority is not None and self.multiplicity != 1:
            self.button["PDR_button"].config(state=tk.NORMAL)
        if PDR_Priority is None:
            self.button["PDR_button"].config(state=tk.DISABLED)
        return (PDR_Priority)

    def update_method_type(self):
        method_type = self.method_type.get()
        if method_type:
            self.method_type_entry = True
            if method_type == 'Equal Bond-Distributed Set':
                self.cheminst_str_type = tk.StringVar(value="Asymmetric")
                self.radio_widgets["Symmetric"].config(state=tk.DISABLED)
            if method_type == 'Rumer Set' or method_type == 'Equal Bond-Distributed Set':
                #self.enable_widgets_in_frame(self.frames["rumer_set_type_frame"])
                #self.disable_widgets_in_frame(self.frames["ChemInst_set_type_frame"])
                #self.disable_widgets_in_frame(self.frames["VB-set_type_frame"])
                self.disable_widgets_in_frame(self.frames["symmetry_set_order_frame"])
                self.disable_widgets_in_frame(self.frames["ovlp_cal_frame"])
                self.disable_widgets_in_frame(self.frames["prio_label_frame"])
                self.disable_widgets_in_frame(self.frames["diff_qualities_frame"])
                self.disable_widgets_in_frame(self.frames["PDB_PDR_Button_frame"])
                self.disable_widgets_in_frame(self.frames["priority_str_frame"])
                #self.disable_widgets_in_frame(self.frames["mout_frame"])
            elif method_type == 'Chemical Insight Set':
                #self.disable_widgets_in_frame(self.frames["rumer_set_type_frame"])
                self.enable_widgets_in_frame(self.frames["ChemInst_set_type_frame"])
                self.enable_widgets_in_frame(self.frames["VB-set_type_frame"])
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
        print('cheminst_type1',cheminst_type)
        if cheminst_type == 'Asymmetric':
            self.button["more_button"].config(state=tk.DISABLED)  # Initially disable the button
            print('cheminst_type2',cheminst_type)
        else:
            self.button["more_button"].config(state=tk.NORMAL) 
            #self.cheminst_str_type_entry = True
            print('cheminst_type3',cheminst_type)
        return (cheminst_type)

    def symmetry_set_order_type_read(self):
        set_order_type = self.symmetry_set_order_type.get()
        if set_order_type:
            self.symmetry_set_order_type_entry = True
            print('set_order_type',set_order_type)
            return (set_order_type)

    def get_mout_number(self):
        mout_entry_window = tk.Toplevel(self.root, padx=10, pady=10)
        mout_entry_window.title("Set Order")
        mout_entry_window.geometry("550x300")
        mout_entry_window.configure(background="lightblue")
        mout_entry_frame = ttk.Frame(mout_entry_window, style="Colour_Frame.TFrame")
        mout_entry_frame.grid(row=0, column=0, sticky=tk.W)
        mout_entry_label= ttk.Label(mout_entry_frame, 
                                            text="If All Sets or Best Sets are opted the number \n"
                                                "of total sets can be millians so users can put a\n"
                                                "maximum number of set (eg. 13 or 1007) or a range \n"
                                                "(eg. 5000-5015). Default maxima sets are 75000 \n"
                                            , style="Colour_Label3.TLabel")
        mout_entry_label.grid(row=0, column=0, columnspan=3, padx=10, pady=10)
        mout_label = ttk.Label(mout_entry_frame, text=" Enter Value: ", style="Colour_Label.TLabel")
        mout_label.grid(row=1, column=0, padx=10, pady=10)
        mout_entry = ttk.Entry(mout_entry_frame, width=10)
        mout_entry.grid(row=1, column=1, padx=10, pady=10)
        enter_button = ttk.Button(mout_entry_frame, text="Enter"
                                  , command=lambda:self.get_maximum_num_output(mout_entry))
        enter_button.grid(row=1, column=3, padx=10, pady=10)
        close_button = ttk.Button(mout_entry_frame, text="Close", command=mout_entry_window.destroy)
        close_button.grid(row=2, column=0, columnspan=3, pady=10)

    def get_set_order(self):
        #if self.get_set_order_window is None:
        self.get_set_order_window = tk.Toplevel(self.root, padx=10, pady=10)
        self.get_set_order_window.title("Set Order")
        self.get_set_order_window.geometry("550x300")
        self.get_set_order_window.configure(background="lightblue")
        symmetry_set_order_frame = ttk.Frame(self.get_set_order_window, style="Colour_Frame.TFrame")
        symmetry_set_order_frame.grid(row=0, column=0, sticky=tk.W)

        symmetry_set_order_label= ttk.Label(symmetry_set_order_frame, 
                                            text="Select any order of arrangement of the symmetric \n"
                                                "group below. It could help searching a symmetric \n"
                                                "set easyer.'Quality' arrange symmetric groups from \n"
                                                "higher quality to lower quality. other two options \n"
                                                "'big to small' and 'small to big' arrange the groups \n"
                                                "according to their sizes. Default is Quality\n"
                                            , style="Colour_Label3.TLabel")
        symmetry_set_order_label.grid(row=0, column=0, columnspan=3, padx=10, pady=10)
        options=("Quality-Arrange", "Big-to-Small", "Small-to-Big")
        for i, option in enumerate(options, start=0):
            button = ttk.Radiobutton(
                symmetry_set_order_frame,
                text=option,
                value=option,
                variable=self.symmetry_set_order_type,
                command=self.symmetry_set_order_type_read,
                style="Custom.TRadiobutton"
            )
            button.grid(row=1, column=i, padx=10, pady=10)
        close_button = ttk.Button(symmetry_set_order_frame, text="DONE", command=self.get_set_order_window.destroy)
        close_button.grid(row=2, column=0, columnspan=3, pady=10)


    def Insert_PDB(self):
        if self.prio_bond_window is None:
            self.prio_bond_window = tk.Toplevel(self.root, padx=10, pady=10)
            self.prio_bond_window.title("Pre-defined Bonds")
            self.prio_bond_window.geometry("450x450")
            self.prio_bond_window.configure(background="lightblue")

            frame = ttk.Frame(self.prio_bond_window, style="Colour_Frame.TFrame")
            frame.grid(row=0, column=0, sticky=tk.W)

            # Scrollable frame for dynamic fields
            canvas = tk.Canvas(self.prio_bond_window, background="lightblue", height=300, width=415)
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

            label = ttk.Label(frame, text="Total Bonds", style="Colour_Label.TLabel")
            label.grid(row=0, column=0, padx=10, pady=10)
            prio_bond_entry = ttk.Entry(frame, width=10)
            prio_bond_entry.grid(row=0, column=1)

            Insert_button = ttk.Button(
                    frame, text="Insert", command=lambda: self.generate_prio_bond_fields
                    (scrollable_frame, prio_bond_entry)
                    )
            Insert_button.grid(row=0, column=2, padx=10, pady=10)
            orb_info_button = ttk.Button(frame, text="Orbital_Info",command=self.Orbital_Info)
            orb_info_button.grid(row=0, column=3, padx=10, pady=10)

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

            prio_bond_entry = ttk.Entry(frame1, width=10)
            prio_bond_entry.grid(row=i, column=1, padx=10, pady=10)

            self.prio_bond_entries.append(prio_bond_entry)

    def get_prio_bond_data(self):
        """
        Retrieve data from the dynamically created structure entry fields.
        """
        self.prio_bonds = []

        for entry in self.prio_bond_entries:
            text = entry.get().strip()
            text = re.sub(r'[-,]', ' ', text)

            if text:
                nums = [int(x) for x in text.split()]
                a, b = nums
                print('nums',nums)
                if len(nums) > 2:
                    tk.messagebox.showerror("Invalid Input",f"Number of orbitals must \n"
                                             "not exceed 2")
                    return
                self.prio_bonds.append(a)
                self.prio_bonds.append(b)
        
        print("Entered bonds:", self.prio_bonds)
        self.prio_bond_window.destroy()
        self.prio_bond_window = None

    def Insert_PDR(self):
        if self.prio_rads_window is None:
            self.prio_rads_window = tk.Toplevel(self.root, padx=10, pady=10)
            self.prio_rads_window.title("Priority Radicals")
            self.prio_rads_window.geometry("650x460")
            self.prio_rads_window.configure(background="lightblue")

            nlp = abs(GlobVar.num_orbital - GlobVar.num_electron)

            frame = ttk.Frame(self.prio_rads_window, style="Colour_Frame.TFrame")
            frame.grid(row=0, column=0, sticky=tk.W)

            # Scrollable frame for dynamic fields
            canvas = tk.Canvas(self.prio_rads_window, background="lightblue", height=350, width=620)
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

           # if nlp == 0:
           #     rad_count = 0
           # else:
            rad_count =1

            Insert_button = ttk.Button(
                    frame, text="Create blank field",
                   # command=lambda: self.generate_prio_rads_fields(scrollable_frame, nlp, rad_count)
                    command=lambda: self.generate_prio_rads_fields(scrollable_frame, rad_count)
                    )
            Insert_button.grid(row=0, column=2, padx=10, pady=10)
            orb_info_button = ttk.Button(frame, text="Orbital_Info",command=self.Orbital_Info)
            orb_info_button.grid(row=0, column=3, padx=10, pady=10)

            close_button = ttk.Button(frame1, text="DONE", command=self.get_prio_rads_data)
            close_button.grid(row=0, column=0, pady=10, columnspan=2)

    def generate_prio_rads_fields(self, frame1, rad_count):
        """
        Generate labels and entry fields dynamically based on user input.
        """
        self.i = self.i + 1
        print('self.i',self.i)
        row_val = self.i
        prio_lnp_label = ttk.Label(
                frame1, text=f"Preferred lone_pair {self.i}:", style="Colour_Label1.TLabel"
                )
        prio_lnp_label.grid(row=row_val, column=0, padx=10, pady=10)

        prio_lnp_entry = ttk.Entry(frame1, width=10)
        prio_lnp_entry.grid(row=row_val, column=1, padx=10, pady=10)
        prio_rads_label = ttk.Label(
                frame1, text=f" and radicals {self.i}:", style="Colour_Label1.TLabel"
                )
        prio_rads_label.grid(row=row_val, column=2, padx=10, pady=10)

        prio_rads_entry = ttk.Entry(frame1, width=20)
        prio_rads_entry.grid(row=row_val, column=3, padx=10, pady=10)

        self.prio_rads_entries.append(prio_rads_entry)
        self.prio_lnp_entries.append(prio_lnp_entry)

    def get_prio_rads_data(self):
        """
        Retrieve data from the dynamically created structure entry fields.
        """
        flg = 0
        for lnp_entry, rad_entry in zip(self.prio_lnp_entries, self.prio_rads_entries):
            val = lnp_entry.get().strip()
            print("val",val)

            if val == "0" and flg == 0:
                text = rad_entry.get().strip()

                if text:
                    # Split by spaces or commas, ignore empty strings
                    for x in re.split(r'[\s,]+', text):
                        if x:
                            self.pref_rad_py.append(int(x))  # convert to int if needed
                    flg = 1

                self.numrad_py = len(self.pref_rad_py)

                print("self.pref_rad_py",self.pref_rad_py)

            else:
            #self.numrad = 0
                self.pref_rad_lnp = []
                numlp = []
                numradlp = []
                self.numudr_py += 1

                lnp_text = lnp_entry.get().strip()
                rad_text = rad_entry.get().strip()

                if lnp_text:
                    self.pref_rad_lnp.extend(lnp_text.split())
                    numlp.append(len(self.pref_rad_lnp))
                    lnp = len(self.pref_rad_lnp)

                self.pref_rad_lnp.append('0')

                if rad_text:
                    self.pref_rad_lnp.extend(rad_text.split())
                    numradlp.append(len(self.pref_rad_lnp)-(1+lnp))

                self.pref_rad_lp_py.append(self.pref_rad_lnp)
                self.numlp_py.append(numlp)
                self.numradlp_py.append(numradlp)
                print("self.pref_rad_lp",self.pref_rad_lp_py)
                print("numrad",self.numrad_py)
                print("numlp_py",self.numlp_py)
                print("numradlp",self.numradlp_py)
        self.prio_rads_window.destroy()
        self.prio_rads_window = None
        #return data
        
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
            orb_info_button = ttk.Button(frame, text="Orbital_Info",command=self.Orbital_Info)
            orb_info_button.grid(row=0, column=3, padx=10, pady=10)

            close_button = ttk.Button(frame1, text="DONE", command=self.get_priority_str_data)
            close_button.grid(row=0, column=0, pady=10, columnspan=2)

    def generate_priority_str_fields(self, frame1, str_entry):
        """
        Generate labels and entry fields dynamically based on user input.
        """
        try:
            self.str_number = int(str_entry.get())
        except ValueError:
            tk.messagebox.showerror("Invalid Input", "Please enter a valid number.")
            return

        # Clear previous fields
        for widget in frame1.winfo_children():
            widget.destroy()

        self.prio_structure_entries.clear()

        for i in range(self.str_number):
            str_label = ttk.Label(frame1, text=f"Structure {i+1}:", style="Colour_Label1.TLabel")
            str_label.grid(row=i, column=0, padx=10, pady=10)

            struc_entry = ttk.Entry(frame1, width=50)
            struc_entry.grid(row=i, column=1, padx=10, pady=10)

            self.prio_structure_entries.append(struc_entry)  # Save reference

    def get_priority_str_data(self):
        """
        Retrieve data from the dynamically created structure entry fields.
        """
        self.strt_strucs = []

        for entry in self.prio_structure_entries:
            text = entry.get().strip()
            text = re.sub(r'[-,]', ' ', text)

            if text:
                nums = [int(x) for x in text.split()]
                if len(nums) > GlobVar.num_orbital:
                    tk.messagebox.showerror("Invalid Input",f"Number of orbitals must \n"
                                             "not exceed {GlobVar.num_orbital}.")
                    return
                self.strt_strucs.append(nums)
        
        print("Entered Structures:", self.strt_strucs)
        self.priority_str_window.destroy()
        self.priority_str_window = None

    def Orbital_Info(self):
        tk.messagebox.showinfo("Orbital Informations",f"there are {GlobVar.num_iao} inactive orbitals\n" 
                                                      f"and {GlobVar.num_orbital} active orbitals,\n" 
                                                      f"which are: {GlobVar.active_orbitals}")


    def get_keywds(self):
        """Retrieve and process default keyword values."""
        def map_string_to_int(value, mapping, default=0):
            """Helper function to map string values to integers."""
            return mapping.get(value, default)

        checksym, nmbond, symm = 0, 0, 1

        method_mapping = {
                'Chemical Insight Set': 1,
                'Rumer Set': 0,
                'Equal Bond-Distributed Set': 2,
                'Check Set if Symmetric': 3
                }
        chinst = map_string_to_int(self.update_method_type(), method_mapping)
        GlobVar.Rum_Ch = chinst

        cheminst_mapping = {'Symmetric': (1, 0), 'Asymmetric': (0, 1)}
        cheminst = self.cheminst_str_type_read()
        symm, asymm = cheminst_mapping.get(cheminst, (0, 1))
        #checksym = 1 if cheminst == 'Checksymm' else 0

        set_order_mapping = {'Quality-Arrange': 0, 'Small-to-Big': 1, 'Big-to-Small': 2}
        set_order = map_string_to_int(self.symmetry_set_order_type_read(), set_order_mapping)

        str_type_mapping = {'Both': 1, 'Covalent': 2, 'Ionic': 3}
        strtype = map_string_to_int(self.Update_Str_Type(), str_type_mapping)
        GlobVar.CovIon = strtype

        settype_mapping = {'Single Set': 0, 'Best Sets': 1, 'All Sets': 2}
        nset = map_string_to_int(self.ChemInst_set_type_read(), settype_mapping)
        GlobVar.num_sets = nset
        if nset == 1:
            nset = 2  # it opt for All Sets as Best Sets will be calculated in the output section

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
        prio_list = [itbp, nnbp, sybp, mnbondp, radicalp]
        count = 0
        for i in range (1,len(prio_list)+1):
            for j in range (len(prio_list)):
                if i == prio_list[j]:
                    count += 1
                    prio_list[j] = count

        itbp, nnbp, sybp, mnbondp, radicalp = prio_list

        #if self.get_PDB_Priority() != 'None':
        #    nmbond = self.bond_number


        if hasattr(self, 'strt_strucs') and hasattr(self, 'str_number'):
            strt_strucs = self.strt_strucs
            str_number = self.str_number
            print("Number of Structures:", str_number)
            print("Structures:", strt_strucs)
        else:
            strt_strucs = 0
            str_number = 0
            
        if hasattr(self, 'prio_bonds') and hasattr(self, 'bond_number'):
            main_bonds = self.prio_bonds
            nmbond = self.bond_number
            print("Number of bonds:", nmbond)
            print("bonds:", main_bonds)
        else:
            main_bonds = 0
            nmbond = 0

        print("I am here 1")
        Run_button.config(state=tk.NORMAL)
        return (chinst, symm, asymm, checksym, strtype, set_order, nset, mout, overlap,
                itbp, nnbp, sybp, mnbondp, radicalp, nmbond, str_number, strt_strucs, main_bonds,
                self.pref_rad_lp_py, self.numlp_py, self.numradlp_py, self.pref_rad_py,
                self.numrad_py, self.numudr_py)
class Help_info:
    def __init__(self):
        pass

    def open_help_window(self):
        help_win = tk.Toplevel()
        help_win.title("Help - Run Options")
        help_win.geometry("900x700")

        # Left menu
        menu_frame = tk.Frame(help_win, width=200, bg='lightgray')
        menu_frame.pack(side=tk.LEFT, fill=tk.Y)

        # Right content area
        content_frame = tk.Frame(help_win)
        content_frame.pack(side=tk.RIGHT, expand=True, fill=tk.BOTH)

        topics = [
            "1. Basic Info",
            "2. Valence Bond Structure",
            "    2.1 Covalent structures",
            "    2.2 Ionic structures",
            "3. Generate structures",
            "  3.1 Basic Informations",
            "  3.2 Orbitals",
            "  3.3 Keywords",
            "4. Runs & Output",
            "5. Visualisation",
            "6. Debug Mode"
        ]

        listbox = tk.Listbox(menu_frame, width=30)
        for topic in topics:
            listbox.insert(tk.END, topic)

        listbox.pack(pady=10, padx=10, fill=tk.BOTH)

        text = tk.Text(content_frame, wrap=tk.WORD)
        text.pack(fill=tk.BOTH, expand=True)

        # Paper URLs
        pub1_url = "https://pubs.acs.org/doi/10.1021/acs.jctc.2c01000"
        pub2_url = "https://pubmed.ncbi.nlm.nih.gov/40511678/"

        # -------------------------
        # BASIC INFO
        # -------------------------

        def show_basic():

            text.config(state=tk.NORMAL)
            text.delete("1.0", tk.END)

            intro = (
            "This software generates, analyzes and visualizes "
            "Heitler-London-Slater-Pauling (HLSP) functions for valence bond calculations"
            "which are also commonly called Valence Bond (VB) structures."
            "This software can provide one or many sets of HLSP functions are with enhanced chemical "
            "interpretation and symmetry-adapted structure sets.\n\n"

            "The methodology implemented in this program is based on "
            "two research publications:\n\n"
            )

            text.insert(tk.END, intro)

            pub1 = (
            "1) New Methodology to Produce Sets of Valence Bond Structures "
            "with Enhanced Chemical Insights,"
            )
            pub1_info1 = " Roy & Shurki, JCTC (2023).\n"

            text.insert(tk.END, pub1, "link1")
            text.insert(tk.END, pub1_info1)

            pub1_info2 = (
            "   • Systematic VB structure generation\n"
            "   • Chemically intuitive VB sets\n\n"
            )

            text.insert(tk.END, pub1_info2)

            text.tag_config("link1", foreground="blue", underline=1)
            text.tag_bind("link1","<Button-1>",lambda e: webbrowser.open(pub1_url))

            pub2 = (
            "2) The Topological Way — A New Methodology to Construct "
            "Symmetric Sets of Valence-Bond Structures,"
            )
            pub2_info1 = (
            " Roy & Shurki, JCP (2025).\n"
            )

            text.insert(tk.END, pub2,"link2", pub2_info1)

            pub2_info2 = (
            "   • Topology-based VB structure generation\n"
            "   • Symmetry-adapted structure sets\n\n"
            )

            text.insert(tk.END,pub2_info2)

            text.tag_config("link2", foreground="blue", underline=1)
            text.tag_bind("link2","<Button-1>",lambda e: webbrowser.open(pub2_url))

            text.config(state=tk.DISABLED)

        # -------------------------
        # VB STRUCTURE
        # -------------------------

        def show_vb():

            text.config(state=tk.NORMAL)
            text.delete("1.0", tk.END)

            vb_text = (

            "Valence Bond (VB) structures represent specific pairing "
            "distributions of localized orbitals.\n\n"

            "A VB structure is expressed as a sequence of integers where "
            "each integer corresponds to an atomic orbital.\n\n"

            "Structure layout:\n\n"

            "   • Inactive orbitals\n"
            "   • Active lone pairs\n"
            "   • Bonded orbital pairs\n"
            "   • Radical orbitals\n\n"

            "Example system:\n"

            "   5 inactive orbitals\n"
            "   1 lone pair\n"
            "   2 singlet bonds\n"
            "   1 radical orbital\n\n"

            "Example structure:\n\n"

            "   1 1 2 2 3 3 4 4 5 5 6 6 7 8 9 10 11\n\n"

            "Compact notation:\n"

            "   1:5 6 6 7 8 9 10 11\n\n"

            "Structure count formula:\n\n"

            "[(2S+1)N! / ((N/2+S+1)! (N/2−S)!)] × [N! / (L!(N−L)!)]\n\n"

            "N : active singly occupied orbitals\n"
            "S : total spin\n"
            "L : number of lone pairs\n"

            "There are two different types of structures for two different bondings"
            "a) Covalent Bonding b) Ionic Bonding"

            )

            text.insert(tk.END,vb_text)
            text.config(state=tk.DISABLED)
        # -------------------------
        # Covalent STRUCTURE
        # -------------------------

        def show_Cov_str():

            text.config(state=tk.NORMAL)
            text.delete("1.0", tk.END)

            vb_text = ("Covalent bonds are singlet pairings of two singly-occupied orbitals, which"
                       "are represented by two consicutive orbital numbers written in the structures."
                        "covalent structures must contain one or more singlet pairs as well as all "
                       "the active orbitals must be includede. As an example, if a system has five"
                       " active orbitals and spin 1/2, then the covalent structure of this system will be\n\n"
                       "                      1 2 3 4 5 or 1 3 2 5 4\n\n"
                       "1 2, 3 4 or 1 3, 2 5 are singlet pairs and 5 and 4 are radicals")

            text.insert(tk.END,vb_text)
            text.config(state=tk.DISABLED)
        # -------------------------
        # IONIC STRUCTURE
        # -------------------------

        def show_Ion_str():

            text.config(state=tk.NORMAL)
            text.delete("1.0", tk.END)

            vb_text = ("Ionic bonds produced when and electron fully shifted towards other atom,"
                       " therefotre one orbital get two electrons or a lone pair and become positively"
                       " charged whereas, other atom posses an empty orbital and becomes a negatively"
                       " charged. Therefore, in the structure we have an orbital with a lone pair as well"
                       " as one orbital is removed for each ionic bond. Now an ionic structure can have"
                       " one ionic bond or can have more. A system must have the identical number of"
                       " covalent and ionic bonds. As an example, if a system has five"
                       " active orbitals and spin 1/2, then the ionic structures of this system will be\n\n"
                       "                         1 1 3 4 5 or 1 1 2 2 4\n\n"
                       "1 1, 2 2 are ionic structure whereas, 3 4 is a singlet pair and 5 and 4 are radicals")

            text.insert(tk.END,vb_text)
            text.config(state=tk.DISABLED)

        # -------------------------
        # str_generation
        # -------------------------

        def str_gen():
            text.config(state=tk.NORMAL)
            text.delete("1.0", tk.END)

            str_gen_text = ("To generate a basic set of structures need to provide a few basic informations"
                       " about the system are listed below:\n\n"
                       " 1: Basic Informations About The Molecular System:\n"
                       "    1.1: Name of the system\n"
                       "    1.2: Geometric informations\n"
                       "         1.2.1: Unit of geometry values\n"
                       "         1.2.2: Brows to upload previously saved geometry file\n"
                       "                or insert geometry manually\n"
                       "    1.3: Active space descriptions\n"
                       "         1.3.1: Number of Active Orbitals (nao)\n"
                       "         1.3.2: Number of Active Electrons (nae)\n"
                       "         1.3.3: Multiplicity of Active Part (nmul)\n"
                       " 2: Orbitals:\n"
                       " 3: Keywords:\n\n")

            text.insert(tk.END,str_gen_text)
            text.config(state=tk.DISABLED)

        # -------------------------
        # BASIC INFO
        # -------------------------

        def basic_info():

            text.config(state=tk.NORMAL)
            text.delete("1.0", tk.END)

            basic_info_text = (

            "\n To generate a set of structures the following system information "
            "must be provided:\n\n"

            "Chemical Formula:\n\n"

            " It must be provided in with there Stoichiometric subscripts.\n\n"

            "Examples:\n"

            "   Benzene → C6H6\n"
            "   Salicylic acid → C7H6O3\n\n"

            "For reactions use reactants or product symbols separated with underscores,\n"
            "if reactants or products are A + B.\n\n"

            "   A_B\n\n"
            "For reaction 2H2 + O2 = 2H2O. It can be written in two ways shown below:\n\n"
            "       1) 2H2_O2  or 2) 2H2O\n\n"

            "\n\nNext we need to provide Geometry informations of the system\n\n"

            "Two options are available to insert Geometries:\n"

            "   1) Upload previously saved geometry file\n"
            "                    or\n"
            "   2) Enter manually\n\n"
            "The geometry file must have a specific format there should be 4 or 5 columns\n"
            "Column two the Atomic Number is optional\n\n"
            "Column 1         Column 2         Column 3         Column 4         Column 5\n"
            "Atomic Formula   Atomic number    X coord          Y coord          Z coord\n\n"
            "After uploding the geometry files it can be seen by clicking 'View_Geometry' button\n"


            "Code supports two different units of geometries as shown below:\n"

            "   • Bohr\n"
            "   • Angstrom\n"
            " User needs to choose one of these option according to there inserted geometry\n\n"

            "\n\nThen information about the active space of the system needs to provude: \n"
            "Active space definition is fully depends upon user's\n"
            "   • Number of Active Orbitals (NAO)\n"
            "   • Number of Active Electrons (NAE)\n"
            "   • Multiplicity (2S+1)\n\n"

            )

            text.insert(tk.END,basic_info_text)
            text.config(state=tk.DISABLED)

        # -------------------------
        # ORBITALS
        # -------------------------

        def show_orb():

            text.config(state=tk.NORMAL)
            text.delete("1.0", tk.END)

            orb_text=(

            "Need to click ORBITALS button after inserting all basic informations\n"
            "about the molecule or reaction. In the orbital panel one needs\n"
            "to provide only active orbital informations. The atom number on which\n"
            "the orbital is situated and the type of the orbitals.\n\n"

            "1) Atom Number\n"
            "2) Orbital Type\n\n"

            "The orbital types can be providede in two different ways: 1) as fragment or\n"
            "hybrid orbital ie, combinations of atomic basis or 2) simply pi1, pi2, pi3 and sigma.\n"
            "They are looks like:\n"

            "   • s or sig (sigma) `\n"
            "   • px or pi1\n"
            "   • py or pi2\n"
            "   • pz or pi3\n"
            )

            text.insert(tk.END,orb_text)
            text.config(state=tk.DISABLED)

        # -------------------------
        # KEYWORDS
        # -------------------------

        def keywd():

            text.config(state=tk.NORMAL)
            text.tag_config("bold", font=("Helvetica", 11, "bold"))
            text.tag_config("underline", underline=1)
            text.delete("1.0", tk.END)

            text.insert(tk.END,"In keyword input pane all the keywords have default values already set up\n"
                        "Therefore one can just press the Insert button and run the program to get\n"
                        "a single symmetric covalent Chemically insightfull set. Bellow all the\n "
                        "Keywords and it's meanings are provided one by one:\n\n")

            text.insert(tk.END,"1.Main Operational Keywords:\n\n","bold")
            text.insert(tk.END,"  1.1","bold") 
            text.insert(tk.END," Chemical Insight Set:",("bold","underline")) 
            text.insert(tk.END," The structures of the set will be selected\n"
                        "      according to the chemical qualities and their priorities. There are\n"
                        "      five different chemical qualities have been considerd among them\n"
                        "      three are natural qualities and two are user preffered qualities\n"
                        "      Natural qualities are explained below one by one:\n\n")
            text.insert(tk.END,"      •","bold") 
            text.insert(tk.END," Intra Atomic Bond (IAB) quality:",("bold","underline")) 
            text.insert(tk.END," it penalize the structures possese\n"
                        "        bonds between two orbitals of the same atom. more IAB more penalty\n\n")
            text.insert(tk.END,"      •","bold") 
            text.insert(tk.END," Near Atomic Bond (NAB) quality:",("bold","underline")) 
            text.insert(tk.END," it penalize the structures possese\n"
                        "        bonds between two atoms situated farther than covalent radius. more\n"
                        "        distant bonds means more penalty. \n\n")
            text.insert(tk.END,"      •","bold") 
            text.insert(tk.END," Symmetry Breaking Bond (SBB) quality:",("bold","underline")) 
            text.insert(tk.END," it penalize the structures possese\n"
                        "        bonds between two orbitals with different symmetries that means bond \n"
                        "        between sigma-sigma or pi1-pi1 is more preferable than bond between \n"
                        "        sigma-pi1 or pi1-pi2. Presence of more this type of bonds get more\n"
                        "        penalty.\n\n")
            text.insert(tk.END,"      •","bold") 
            text.insert(tk.END," Pre-Defined Bond (PDB) quality:",("bold","underline")) 
            text.insert(tk.END," it favour the structures possese\n"
                        "        bonds pre-defined by the user. Presence of more this type of bonds in\n"
                        "        the structures get more favour.\n\n")
            text.insert(tk.END,"      •","bold") 
            text.insert(tk.END," Pre-Defined Radical (PDR) quality:",("bold","underline")) 
            text.insert(tk.END," it favour the structures possese\n"
                        "        radicals pre-defined by the user. Presence of more this type of\n"
                        "        radicals in the structures get more favour.\n\n\n")

            text.insert(tk.END,"  1.2","bold") 
            text.insert(tk.END," Rumer Set:",("bold","underline")) 
            text.insert(tk.END," The structures in the set are selected according to Rumer's\n"
                        "      methodology. In this approach, a linearly independent set of valence\n"
                        "      bond structuresis constructed using a diagrammatic representation.\n"
                        "      First, all singly occupied active orbitals are arranged on a circle.\n"
                        "      Bonds are then drawn between pairs of orbitals.\n"
                        "      A valid linearly independent structure is obtained when the bonds\n"
                        "      do not cross each other.\n\n")

            text.insert(tk.END,"  1.3","bold") 
            text.insert(tk.END," Equal Bond-Distributed Set:",("bold","underline")) 
            text.insert(tk.END," In a whole set of structures many bonds can be\n"
                        "      repeated multiple times. This option will provide sets that contains\n"
                        "      all bonds repeated equal number of times.\n"
                        "      This type of sets can have some important use.\n\n")

            text.insert(tk.END,"  1.4","bold") 
            text.insert(tk.END," Check Set if Symmetric:",("bold","underline")) 
            text.insert(tk.END," This option provide the information about a set\n"
                        "      if it is symmetric or not. This option needs sets to be uploaded.\n\n")

            text.insert(tk.END,"2.General Keywords:\n\n","bold")
            text.insert(tk.END,"  2.1","bold") 
            text.insert(tk.END," HLSP Function Set Type:",("bold","underline")) 
            text.insert(tk.END," Two different types of HLSP functions set can be \n"
                        "      generated which are, symmetric and asymmetric. Here we are only \n"
                        "      considering symmetry of the active part of the system. If one select\n"
                        "      Symmetric, then they can choose one more option which is basically to\n"
                        "      to choose the order of the symmetric HLSP function groups. there are\n"
                        "      one of three orders can be chosen and according to that we have more\n"
                        "      and more chances to get that type of structures in the set. The options\n"
                        "      are explained below:\n")
            text.insert(tk.END,"    a)","bold") 
            text.insert(tk.END," Quality-Arrange:",("bold","underline")) 
            text.insert(tk.END," Every HLSP function in a symmetric group\n" 
                        "      has the same Chemical qualities and that represents the chemical quality\n"
                        "      of that entire group. In this option program will try to provide a set\n"
                        "      with more higher chemical quality structures.\n")
            text.insert(tk.END,"    b)","bold") 
            text.insert(tk.END," Big-to-Small:",("bold","underline")) 
            text.insert(tk.END," The symmetric groups can have a single\n" 
                        "      HLSP function or can have multiple. In this option program will try to\n"
                        "      provide a set contains as much as possible bigger or more populated\n"
                        "      groups.")
            text.insert(tk.END,"    c)","bold") 
            text.insert(tk.END," Small-to-Big:",("bold","underline")) 
            text.insert(tk.END," The symmetric groups can have a single\n" 
                        "      HLSP function or can have multiple. In this option program will try to\n"
                        "      provide a set contains as much as possible smaller or less populated groups.\n\n")
            text.insert(tk.END,"  2.2","bold") 
            text.insert(tk.END," HLSP Function Type:",("bold","underline")) 
            text.insert(tk.END," Two types of HLSP functions represents two\n"
                        "      different types of bonding which are, Covalent and Ionic. One of \n"
                        "      either or both can be needed depending on the methodology of\n"
                        "      calculation. One need to select any one of the options among Covalent\n"
                        "      , Ionic or Both.\n")
            text.insert(tk.END,"  2.3","bold") 
            text.insert(tk.END," Set Type:",("bold","underline")) 
            text.insert(tk.END," The possible number of HLSP functions\n"
                        "      is more than the number required to create an independent set.  \n"
                        "      This program can generate all possible number of sets. one can asked\n"
                        "      for a single best set (best according to the chemical qualities) or \n"
                        "      all best sets or all sets or any number of sets. If they want few \n"
                        "      number of sets, they must select 'All Sets' option and then click \n"
                        "      more option and put the number how many sets they want. The different\n"
                        "      options are shown below\n\n")
            text.insert(tk.END,"    a)","bold") 
            text.insert(tk.END," Single Set:",("bold","underline")) 
            text.insert(tk.END," Select this option to get one best\n" 
                        "      set according to the chemical qualities and their priorities are \n"
                        "      given.\n")
            text.insert(tk.END,"    b)","bold") 
            text.insert(tk.END," All Sets:",("bold","underline")) 
            text.insert(tk.END," Select this option to get all possible\n" 
                        "      sets, where the first set will be the best set. This all possible \n"
                        "      set will contain all symmetric and asymmetric set. One thing must\n"
                        "      need to keep in mind that as the number of active orbitals of the \n"
                        "      system will increase the number of set will increase exponentially\n"
                        "      So it is allways safe to provide a realistic number inside the option\n"
                        "      'more'.\n")
            text.insert(tk.END,"    a)","bold") 
            text.insert(tk.END," Best Sets:",("bold","underline")) 
            text.insert(tk.END," In a specific chemical quality and priority\n"
                        "      the system can have more than one best set, and selecting this \n" 
                        "      option one would get all best quality sets. \n\n")

            text.insert(tk.END,"3.Priorities of Chemical Qualities:","bold")
            text.insert(tk.END," There are five chemical qualities what \n"
                        "      user can prioritize these chemical qualities. All five qualities\n "
                        "      may not be available for that particular active space or system\n"
                        "      such as two atomic molecules do not have NAB quality. Two qualities\n"
                        "      can have same priorities. The default priorities are 1st:IAB, 2nd:NAB,\n"
                        "      3rd:SBB, PDB and PDR set in None. To put priority for PDB and PDR one\n"
                        "      need to insert there prefered bonds and prefered radicles through\n"
                        "      insert buttons 'Insert PDB' & 'Insert PDR'\n")
            text.insert(tk.END,"   3.1","bold") 
            text.insert(tk.END," Insert PDB:",("bold","underline")) 
            text.insert(tk.END," To insert PDBs in the opend window 'Pre-Defined Bonds', first need \n"
                        "      to insert total numberof bonds and then put separately in each entry \n"
                        "      boxes. One can insert bonded orbitals separated with coma, dash or \n"
                        "      blank pace and finally click the insert button")
            text.insert(tk.END,"   3.2","bold") 
            text.insert(tk.END," Insert PDR:",("bold","underline")) 
            text.insert(tk.END," To insert PDRs in the opend window 'Pre-Defined Radicals', first need \n"
                        "      to click the button 'Create blank field'. Every click will open two \n"
                        "      entry boxes in a row. First box is for lone pairs and t0he second one \n"
                        "      is for radicals. Two types of situation can appeared for two types of \n"
                        "      HLSP function sets\n" 
                        "        1) with lone pairs \n"
                        "        2) without lone pairs"
                        "      irrespective of covalent or ionic functions. If it is a set without lone\n"
                        "      pairs, one must need to put 0 (zero) in the box for 'preffered lone_pair n'.\n"
                        "      Radicals for without lone pair set can be provided onece only other wise\n"
                        "      only first one will be taken. For sets with lone-pairs one can put as\n"
                        "      many lone pairs as they want. If some lone pair sets are not been given any\n"
                        "      radical prefference the set will be guided by rest of the qualities and if\n"
                        "      no other qualities are selected then the choice of that set would be fully\n"
                        "      arbitrary.")
            text.insert(tk.END,"4. Preferred Structures:","bold")
            text.insert(tk.END," User can insert their preferred \n"
                        "      structures using 'Insert Preferred Structure' button. Programe will try to\n "
                        "      include these structures without judging their qualities (they might or might\n"
                        "      not mached with the preferred qualities of the run). However, it cannot be\n"
                        "      guaranteed the inclusion of all preferred structures as long all of them\n"
                        "      are mutually linearly independent.")
            #text.insert(tk.END,keywd_text)
            text.config(state=tk.DISABLED)

        # -------------------------
        # RUN
        # -------------------------

        def show_run():

            text.config(state=tk.NORMAL)
            text.delete("1.0", tk.END)
            text.insert(tk.END,"  4.1 Various Runs:","bold")
            text.insert(tk.END," The programe can capable \n"
                        "      to produce Chemical Insight Sets, Rumer sets, Equal Bond Distributed Sets\n"
                        "      with various manupulations.")

            text.insert=(

            "Press the RUN button to start the calculation.\n\n"

            "The program generates:\n"

            "   • VB structures\n"
            "   • symmetry adapted sets\n"
            "   • analysis output\n"

            )

            #text.insert(tk.END,run_text)
            text.config(state=tk.DISABLED)

        # -------------------------
        # VISUALISATION
        # -------------------------

        def show_vis():

            text.config(state=tk.NORMAL)
            text.delete("1.0", tk.END)

            vis_text=(

            "Visualization module displays:\n\n"

            "   • molecular orbitals\n"
            "   • electron density\n"
            "   • geometry structures\n"

            )

            text.insert(tk.END,vis_text)
            text.config(state=tk.DISABLED)

        # -------------------------
        # DEBUG
        # -------------------------

        def show_debug():

            text.config(state=tk.NORMAL)
            text.delete("1.0", tk.END)

            debug_text=(

            "Debug mode prints intermediate variables\n"
            "and diagnostic information.\n\n"

            "Useful for:\n"

            "   • debugging inputs\n"
            "   • checking algorithm steps\n"
            "   • development testing\n"

            )

            text.insert(tk.END,debug_text)
            text.config(state=tk.DISABLED)

        # -------------------------
        # MENU HANDLER
        # -------------------------

        def on_select(event):

            selection = listbox.curselection()
            if not selection:
                return

            topic = listbox.get(selection[0])

            if topic=="1. Basic Info":
                show_basic()

            elif topic=="2. Valence Bond Structure":
                show_vb()

            elif topic=="    2.1 Covalent structures":
                show_Cov_str()

            elif topic=="    2.2 Ionic structures":
                show_Ion_str()

            elif topic=="3. Generate structures":
                str_gen()

            elif topic=="  3.1 Basic Informations":
                basic_info()

            elif topic=="  3.2 Orbitals":
                show_orb()

            elif topic=="  3.3 Keywords":
                keywd()

            elif topic=="4. Run & Output":
                show_run()

            elif topic=="5. Visualisation":
                show_vis()

            elif topic=="6. Debug Mode":
                show_debug()

        listbox.bind("<<ListboxSelect>>",on_select)

class drawing_molecule:
    def __init__(self):
        self.atoset = GlobVar.Gl_atoset.copy()
        self.active_atoms = GlobVar.Gl_activeatoms.copy()
        self.atoms = GlobVar.all_atoms.copy()
        self.local_orbital_data = []
        self.local_orbital_data = GlobVar.orbital_data.copy()
        self.nlp_lst = []
        self.rad_lst = []
        self.bond_lst = []
        self.str_number = None
        self.set_number = None
        print('self .active_atoms',self.active_atoms)
        print('self.atoms',self.atoset)

    def str_info(self, nlp_lst, rad_lst, bond_lst, str_number, set_number):
        self.nlp_lst = nlp_lst
        self.rad_lst = rad_lst
        self.bond_lst = bond_lst
        self.str_number = str_number
        self.set_number = set_number

    def VB_view(self):
        if GlobVar.geometry_inserted is True:
            mol_display = tk.Toplevel()
            if GlobVar.molecule_string and not self.str_number:
                mol_display.title(f"{GlobVar.molecule_string}")
            elif GlobVar.molecule_string and self.str_number:
                mol_display.title(f"Molecule: {GlobVar.molecule_string}, Set_number: {self.set_number}, Structure Number: {self.str_number}")
            else:
                mol_display.title("-Molecule-")

            mol_display.geometry("600x700")
            mol_display.configure(background="lightblue")

            mol_frame = tk.Frame(mol_display)
            mol_frame.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

            button_frame = tk.Frame(mol_frame)
            button_frame.pack(fill=tk.BOTH, expand=True)

            center_frame = tk.Frame(button_frame)
            center_frame.pack(expand=True)

            canvas_frame = tk.Frame(mol_frame)
            canvas_frame.pack(fill=tk.BOTH, expand=True)

            box_frame = tk.Frame(mol_frame, width=100, height=100, bg='blue')
            box_frame.pack(pady=10)
            box_frame.pack_propagate(False)
            
            fig = plt.Figure(figsize=(6, 6))
            fig.patch.set_facecolor('lightblue')
            ax = fig.add_subplot(111, projection='3d')
            ax.set_facecolor('lightblue')
            
            # drawing axis
            ax.plot([0, 2], [0, 0], [0, 0], color='green', linewidth=2, linestyle=':', alpha=0.8)
            ax.plot([0, 0], [0, 2], [0, 0], color='green', linewidth=2, linestyle=':', alpha=0.8)
            ax.plot([0, 0], [0, 0], [0, 2], color='green', linewidth=2, linestyle=':', alpha=0.8)
            ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
            # Add axis labels
            ax.text(2, 0, 0, 'X', color='green', fontsize=14, fontweight='bold', alpha=0.8)
            ax.text(0, 2, 0, 'Y', color='green', fontsize=14, fontweight='bold', alpha=0.8)
            ax.text(0, 0, 2, 'Z', color='green', fontsize=14, fontweight='bold', alpha=0.8)

            # drawing atoms
            old_colors = []
            self.selected_atoms = []
            atom_scatters = []
            scale_factor = 0.5
            for entry in self.atoms:
                if entry['sl_num'] in self.active_atoms:
                    x1 = entry['x']
                    y1 = entry['y']
                    z1 = entry['z']
                    symbol = entry['atom']

                    radius = GlobVar.at_covrad.get(symbol)  # fallback value if element missing
                    size = (radius * scale_factor)**2
                    print('size of the atom',size)

                    scatters = ax.scatter(x1, y1, z1,                   # drawing atoms
                               color=GlobVar.colors.get(symbol, 'green'),
                               s=size, alpha=0.4)

                    atom_scatters.append(scatters)
                    ax.text(x1 - 0.15, y1 - 0.15, z1 - 0.15, symbol, 
                            color='black', fontsize=10, ha='center', va='bottom') # puttong names
                    clr = GlobVar.colors.get(symbol)
                    old_colors.append(clr)
            #print('scatters',scatters)

            # drawing bonds
            for entry1, entry2 in combinations(self.atoms, 2):
                if entry1['sl_num'] in self.active_atoms and entry2['sl_num'] in self.active_atoms:
                    x1, y1, z1 = entry1['x'], entry1['y'], entry1['z']
                    x2, y2, z2 = entry2['x'], entry2['y'], entry2['z']

                    distance = math.sqrt(
                            (x2-x1)**2 +
                            (y2-y1)**2 +
                            (z2-z1)**2
                            )

                    r1 = GlobVar.at_covrad.get(entry1['atom'])
                    r2 = GlobVar.at_covrad.get(entry2['atom'])

                    if r1 is not None and r2 is not None:
                        covrad = (r1 + r2)/100.0
                        if distance <= covrad: # cheking if the distances between two atoms are bonded
                            ax.plot([x1, x2], [y1, y2], [z1, z2], color='black', linestyle=':', linewidth=2.5, alpha=0.5)
                    else:
                        messagebox.showerror("Not Found",f"covalent radious of {atom1} and/or {atom2}\n"
                                             "are not found in the list")
                            
            
            # Make panes transparent (no background)
            ax.xaxis.set_pane_color((1, 1, 1, 0))
            ax.yaxis.set_pane_color((1, 1, 1, 0))
            ax.zaxis.set_pane_color((1, 1, 1, 0))
            
            # Hide grid lines
            ax.grid(False)
            
            # Hide tick marks and labels
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            
            # Hide box edges (axis lines)
            try:
                ax.w_xaxis.line.set_color((0, 0, 0, 0))
                ax.w_yaxis.line.set_color((0, 0, 0, 0))
                ax.w_zaxis.line.set_color((0, 0, 0, 0))
            except AttributeError:
                # For newer Matplotlib versions
                for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
                    axis.line.set_alpha(0)

            orbital_coord = []
            single_orbs = ()
#            local_orbital_data = []
            print('GlobVar_orbital_data',GlobVar.orbital_data)
            print('local_orbital_data:', self.local_orbital_data)

            # create a list of orbital number and its coordinates to
            # which it is associated
            for i, row in enumerate(self.atoset):
                nonzero_indices = np.nonzero(row)[0]
                nonzero_values = row[nonzero_indices]
                if nonzero_values.size > 0:
                    for entry1 in nonzero_values:
                        print('nonzero_value', nonzero_values)
                        for entry2 in self.atoms:
                            if entry2['sl_num'] == i+1:
                                x = entry2['x']
                                y = entry2['y']
                                z = entry2['z']
                                symbol = entry2['atom']
                                print('atom:',symbol, x, y, z)
                                print('local_orb',self.local_orbital_data)
                        for entry3 in self.local_orbital_data:
                            print('local:',self.local_orbital_data)
                            if entry3['atom_number'] == i+1:
                                orb_type = entry3['orbital_type']
                                print('orb_type',orb_type)

                                single_orbs = (symbol, entry1, orb_type, x, y, z)
                                print('single_orbs',single_orbs)
                                orbital_coord.append(single_orbs)
                                self.local_orbital_data.remove(entry3)
                                #local_orbital_data.remove(entry3)
                                break
                        continue
            print('orbital_data', GlobVar.orbital_data)
            print ('orbital_coord',orbital_coord)

            # print(GlobVar.orbital_data)
            # print(orbital_coord)
            vb_bonds_points = []
            guiding_points = []
            for entry in orbital_coord:
                symbol, orb_num, orb_typ, x0, y0, z0 = entry
                if 's' == orb_typ.lower() or 'sig' in orb_typ.lower():
                    radius = GlobVar.at_covrad.get(symbol)  # fallback value if element missing
                    n = radius/280
                    a, b = n, n + 0.07  # axes lengths of lower part
                    u = np.linspace(0, 2 * np.pi, 40)
                    v = np.linspace(0, np.pi, 40)  # upper hemisphere
                    u, v = np.meshgrid(u, v)
                    
                    # Parametric surface
                    x = a * np.cos(u) * np.sin(v)+x0
                    y = a * np.sin(u) * np.sin(v)+y0
                    z = b * np.cos(v) + z0
                    centers = (x0, y0, z0)
                    centers1 = (orb_num, x0, y0, z0)
                    guiding_points.append(centers)
                    vb_bonds_points.append(centers1)

                    ax.plot_surface(x, y, z, color='red', alpha=0.5) # transparency alpha
                    ax.text(x0, y0 + 0.1, z0 + b + 0.1, orb_num, color='black', fontsize=10, ha='center', va='bottom')

                else:
                    x_u, y_u, z_u, x_l, y_l, z_l = self.orbital_builder(
                            symbol, orb_num, orb_typ, x0, y0, z0)
                    xu_c, yu_c, zu_c = np.mean(x_u), np.mean(y_u), np.mean(z_u) # center of
                    xl_c, yl_c, zl_c = np.mean(x_l), np.mean(y_l), np.mean(z_l) # the surface
                    centers_u = (xu_c, yu_c, zu_c)
                    centers_l = (xl_c, yl_c, zl_c)
                    guiding_points.append(centers_u)
                    guiding_points.append(centers_l)

                    ax.plot_surface(x_u, y_u, z_u, color='red', alpha=0.5) # transparency alpha
                    ax.plot_surface(x_l, y_l, z_l, color='red', alpha=0.5) # transparency alpha

                    if 's' in orb_typ.lower():
                        i, j = x_l.shape[0] // 2, x_l.shape[1] // 2
                        p = np.array([x_l[i, j], y_l[i, j], z_l[i, j]])
                        centers1 = (orb_num, xl_c, yl_c, zl_c)
                        if 'pz' in orb_typ.lower() and '-pz' not in orb_typ.lower():
                            i, j = x_u.shape[0] // 2, x_u.shape[1] // 2
                            p = np.array([x_u[i, j], y_u[i, j], z_u[i, j]])
                            centers1 = (orb_num, xu_c, yu_c, zu_c)
                        vb_bonds_points.append(centers1)
                    else:
                        i, j = x_u.shape[0] // 2, x_u.shape[1] // 2
                        p = np.array([x_u[i, j], y_u[i, j], z_u[i, j]])
                        centers1 = (orb_num, xu_c, yu_c, zu_c)
                        vb_bonds_points.append(centers1)

                    
                    # Compute radial direction from center
                    center = np.array([x0, y0, z0])
                    normal = p - center
                    normal = normal / np.linalg.norm(normal)
                    
                    # Offset the point in front of the surface
                    d = 0.2  # small distance in front
                    p_front = p + d * normal

                    ax.text(
                        p_front[0], p_front[1], p_front[2],  # x, y, z coordinates
                        str(orb_num),                        # Convert orb_num to string (if it's numeric)
                        color='black',
                        fontsize=10,
                        ha='center',
                        va='bottom'
                    )

            # drwaing vb bonds after generating structures
            if len(self.bond_lst) != 0:
                for a, b in self.bond_lst:
                    print('a, b', a, b)
                    for orb_num1, x1, y1, z1 in vb_bonds_points:

                            for orb_num2, x2, y2, z2 in vb_bonds_points:
                                    if int(a) == int(orb_num1) and int(b) == int(orb_num2):
                                        p1 = (x1, y1, z1)
                                        p2 = (x2, y2, z2)
                                        filtered_points = guiding_points.copy()
                                        if p1 in filtered_points:
                                            filtered_points.remove(p1)
                                        if p2 in filtered_points:
                                            filtered_points.remove(p2)

                                        #print('guiding_points',guiding_points)
                                        #print('points',p1, p2)
                                        #print('filtered_points',filtered_points)
                                        line_pts = self.get_line_points(p1, p2, 50)
                                        clear, min_dist = self.is_arc_clear(line_pts, filtered_points, d_min=0.01)
                                        print('clear',clear, min_dist)
                                        if clear:
                                            ax.plot([x1, x2], [y1, y2], [z1, z2], color='blue', linewidth=2.5)
                                        else:
                                            height = 0.1
                                            normal = [1, 1, 1]
                                            #success = False
                                            for i in range(20):
                                                arc_pts = self.generate_arc(p1, p2, normal, height)
                                                clear, min_dist = self.is_arc_clear(arc_pts, filtered_points,
                                                                                    d_min=0.01)
                                                if clear:
                                                    #success = True
                                                    ax.plot(arc_pts[:,0], arc_pts[:,1], arc_pts[:,2], color='blue', linewidth=2.5)
                                                    break
                                                height += 0.1

            if len(self.rad_lst) != 0:
                for a in self.rad_lst:
                    for orb_num, x1, y1, z1 in vb_bonds_points:
                        if int(a) == int(orb_num):
                            ax.scatter(x1, y1, z1, color='blue', s=50)


            if len(self.nlp_lst) != 0:
                for a in self.nlp_lst:
                    for orb_num, x1, y1, z1 in vb_bonds_points:
                        if int(a) == int(orb_num):
                            ax.scatter(x1, y1, z1, color='blue', s=50)
                            ax.scatter(x1, y1, z1+0.5, color='blue', s=50)

            canvas = FigureCanvasTkAgg(fig, master=canvas_frame)
            canvas_widget = canvas.get_tk_widget()
            canvas_widget.pack(fill=tk.BOTH, expand=True)

            close_btn = tk.Button(box_frame, text="Close", command=mol_display.destroy)
            close_btn.pack(expand=True)

    def orbital_builder(self, symbol, orb_num, orb_typ, x0, y0, z0 ):
        print('------',symbol, orb_num, orb_typ,'------')
        print('extended_orbital_coord',x0, y0, z0)
        sflg = 0
        x1, y1, z1 = x0, y0, z0
        # create one ellipsoid for pi orbitals and sphere for 
        # sigma orbitals asssociated with each atoms
        # to make it more looks like an pi orbital add two hulf
        # upper part in spherical and lower part is ellipsoid

        # creating axis points at 0.4 points away from the centre of the atom 
        if 's' in orb_typ.lower():
            print('entered in S... orbs')
            radius = GlobVar.at_covrad.get(symbol)  
            n = radius/250

            a, b = n, 2*n  # axes lengths of lower part

            u = np.linspace(0, 2 * np.pi, 40)
            # making upper part
            v_upper = np.linspace(0, np.pi/2, 20)
            x_upper = a * np.outer(np.cos(u), np.sin(v_upper))
            y_upper = a * np.outer(np.sin(u), np.sin(v_upper))
            z_upper = a * np.outer(np.ones_like(u), np.cos(v_upper))

            ## making lower part
            v_lower = np.linspace(np.pi/2, np.pi, 20)
            x_lower = a * np.outer(np.cos(u), np.sin(v_lower))
            y_lower = a * np.outer(np.sin(u), np.sin(v_lower))
            z_lower = b * np.outer(np.ones_like(u), np.cos(v_lower))

            sflg = 1
            pzflg = 0
            pxflg = 0
            pyflg = 0
            if '-pz' in orb_typ.lower():
                pzflg = 1
                return (x_upper+x0, y_upper+y0, z_upper+z0, x_lower+x0, y_lower+y0, z_lower+z0)
            if 'pz' in orb_typ.lower() and pzflg == 0:
                z_upper = b * np.outer(np.ones_like(u), np.cos(v_upper))
                z_lower = a * np.outer(np.ones_like(u), np.cos(v_lower))
                return (x_upper+x0, y_upper+y0, z_upper+z0, x_lower+x0, y_lower+y0, z_lower+z0)
            if '-px' in orb_typ.lower():
                print('here3')
                x1 = x1 + 0.4
                pxflg = 1
            if 'px' in orb_typ.lower() and pxflg == 0:
                print('here4')
                x1 = x1 - 0.4
            if '-py' in orb_typ.lower():
                print('here5')
                y1 = y1 + 0.4
                pyflg = 1
            if 'py' in orb_typ.lower() and pyflg ==0:
                print('here6')
                y1 = y1 - 0.4
        else:
            print('entered in non  S orbs')
            a, b = 0.2, 0.4  # axes lengths of lower part

            radius = GlobVar.at_covrad.get(symbol)  
            size = radius/100

            u = np.linspace(0, 2 * np.pi, 40)

            # making upper part
            v_upper = np.linspace(0, np.pi/2, 20)
            x_upper = a * np.outer(np.cos(u), np.sin(v_upper))#+x0
            y_upper = a * np.outer(np.sin(u), np.sin(v_upper))#+y0
            z_upper = a * np.outer(np.ones_like(u), np.cos(v_upper))#+z0#+size/1.6

            # making lower part
            v_lower = np.linspace(np.pi/2, np.pi, 20)
            x_lower = a * np.outer(np.cos(u), np.sin(v_lower))#+x0
            y_lower = a * np.outer(np.sin(u), np.sin(v_lower))#+y0
            z_lower = b * np.outer(np.ones_like(u), np.cos(v_lower))#+z0#+size/1.6

            pxflg = 0
            pyflg = 0
            pzflg = 0
            if '-pz' in orb_typ.lower():
                z1 = z1 - 0.4
                pzflg = 1
            if 'pz' in orb_typ.lower() and pzflg == 0:
                z1 = z1 + 0.4
            if '-px' in orb_typ.lower():
                x1 = x1 - 0.4
                pxflg = 1
            if 'px' in orb_typ.lower() and pxflg == 0:
                x1 = x1 + 0.4
            if '-py' in orb_typ.lower():
                y1 = y1 - 0.4
                pyflg = 1
            if 'py' in orb_typ.lower() and pyflg ==0:
                y1 = y1 + 0.4

        # find the direction of the preffered axis 
        p0 = np.array([x0, y0, z0])
        p1 = np.array([x1, y1, z1])
        direction = p1 - p0        #direction = (x1−x0,y1−y0,z1−z0)
        print('direction', direction, p0, p1)
        b = direction / np.linalg.norm(direction)  # normalize to get direction unit-vector
        print('unit_vector', b, np.linalg.norm(direction))

        """ Find the rotation matrix that aligns vec1 to the direction of the preffered axis """
        vec1 = np.array([0, 0, 1]) # original directtion of the orbital
        a = (vec1 / np.linalg.norm(vec1)).reshape(3) # unit vector of vec1
        v = np.cross(a, b) #The cross product gives the axis about which to rotate to align a with b
        c = np.dot(a, b) #cosine of the angle between the two vectors.

        # If the vectors are already aligned, the rotation is the identity matrix (no rotation needed).
        if np.isclose(c, 1):
            R = np.eye(3) # 3X3 Identity matrix

        #If the vectors are exactly opposite, rotate 180° about any axis perpendicular to a.
        elif np.isclose(c, -1):
            orth = np.array([1, 0, 0]) if not np.isclose(a[0], 1) else np.array([0, 1, 0])
            v = np.cross(a, orth)
            v = v / np.linalg.norm(v)
            H = np.array([[0, -v[2], v[1]],[v[2], 0, -v[0]],[-v[1], v[0], 0]])
            R = -np.eye(3) + 2 * np.outer(v, v)

        else:
            # For all other cases, use Rodrigues’ rotation formula to construct the rotation matrix.
            s = np.linalg.norm(v) # sin of the angle between two vectors
            kmat = np.array([[0, -v[2], v[1]],
                         [v[2], 0, -v[0]],
                         [-v[1], v[0], 0]])
            R = np.eye(3) + kmat + kmat @ kmat * ((1 - c) / (s ** 2)) # Rodrigues Formula

        x_upper_r, y_upper_r, z_upper_r = self.rotate_and_translate(x_upper, y_upper, z_upper, R, p0)
        x_lower_r, y_lower_r, z_lower_r = self.rotate_and_translate(x_lower, y_lower, z_lower, R, p0)

        if sflg == 0:
            x_upper_r, y_upper_r, z_upper_r, x_lower_r, y_lower_r, z_lower_r= self.position_scaling(
                x_upper_r, y_upper_r, z_upper_r, x_lower_r, 
                y_lower_r, z_lower_r, x0, y0, z0, orb_typ,  size, 'orb')
        return (x_upper_r, y_upper_r, z_upper_r, x_lower_r, y_lower_r, z_lower_r)

    def position_scaling(self, x_u, y_u, z_u, x_l, 
                         y_l, z_l, x0, y0, z0, orb_typ, size, flg):
        # scaling for shifting the orbitals outside of the surface
        a = np.array([x0, y0, z0])
        x_c = np.mean(x_l)
        y_c = np.mean(y_l)
        z_c = np.mean(z_l)
        b = np.array([x_c, y_c, z_c])
        c = np.array([x_u, y_u, z_u])

        direction = b - a
        norm = np.linalg.norm(direction)

        if flg == 'text':
            d = norm - size*1.15 + 0.1  # distance from the center of atom + an offset for the numbers
        elif flg == 'orb':
            d = norm - size*1.15 # distance from the center of atom
        else:
            d = 0.0
        
        if norm == 0:
            raise ValueError("Zero-length direction vector. Check position scaling function")

        unit_vector = direction / norm

        # Translate point b to be at distance d from a
        offset =  d * unit_vector
        x_upper_r = x_u + offset[0]
        y_upper_r = y_u + offset[1]
        z_upper_r = z_u + offset[2]
        x_lower_r = x_l + offset[0]
        y_lower_r = y_l + offset[1]
        z_lower_r = z_l + offset[2]
        return (x_upper_r, y_upper_r, z_upper_r, x_lower_r, y_lower_r, z_lower_r)

    def rotate_and_translate(self, x, y, z, R, center):
        pts = np.stack([x.flatten(), y.flatten(), z.flatten()])
        pts_rot = R @ pts
        pts_rot[0] += center[0]
        pts_rot[1] += center[1]
        pts_rot[2] += center[2]
        return (pts_rot[0].reshape(x.shape), pts_rot[1].reshape(y.shape), 
                pts_rot[2].reshape(z.shape))

    def generate_arc(self, P1, P2, normal, height=1.0, n_points=100):
        P1, P2, normal = map(np.array, (P1, P2, normal))
        v = P2 - P1
        v_hat = v / np.linalg.norm(v)
        print('v:',v, 'v_hat:',v_hat)
        w = np.cross(normal, v_hat)
        print('wwww',w)
        w_hat = w / np.linalg.norm(w)
        print('W_hat',w_hat)

        M = (P1 + P2) / 2
        C = M + height * w_hat  # center of the arc

        r = np.linalg.norm(P1 - C)
        theta1 = np.arctan2(np.dot(P1 - C, w_hat), np.dot(P1 - C, v_hat))
        theta2 = np.arctan2(np.dot(P2 - C, w_hat), np.dot(P2 - C, v_hat))
        thetas = np.linspace(theta1, theta2, n_points)

        arc_points = []
        for theta in thetas:
            point = C + r * (np.cos(theta) * v_hat + np.sin(theta) * w_hat)
            arc_points.append(point)
        return np.array(arc_points)

    def is_arc_clear(self, arc_points, stored_points, d_min):
        from scipy.spatial.distance import cdist
        distances = cdist(arc_points, stored_points)
        min_dist = np.min(distances)
        return min_dist > d_min, min_dist

    def get_line_points(self, p1, p2, num_points=100):
        """Generate num_points along a line from p1 to p2."""
        p1 = np.array(p1)
        p2 = np.array(p2)
        t = np.linspace(0, 1, num_points)
        points = (1 - t[:, None]) * p1 + t[:, None] * p2  # shape: (num_points, 3)
        return points


class Run_Fort:
    def __init__(self, root, ctrl_class):
        self.root = root
        self.ctrl_class = ctrl_class

    def get_keywds(self, keywd_class):
        self.chinst, self.symm, self.asymm, self.checksym, self.strtype, self.set_order, \
        self.nset, self.mout, self.ovlp, self.itb, self.nnb, self.syb, \
        self.mnbond, self.radical, self.nmbond, self.nstrt1, self.strt_struc1,\
        self.main_bond, self.pref_rad_lp, self.numlp, self.numradlp,\
        self.pref_rad, self.numrad_py, self.numudr_py = keywd_class.get_keywds()

        GlobVar.symm_key = self.symm

        self.strt_struc = np.array(self.strt_struc1, dtype=np.int32, order='F')
        self.main_bond_py = np.array(self.main_bond, dtype=np.int32, order='F')
        self.pref_rad_lp_py = np.array(self.pref_rad_lp, dtype=np.int32, order='F')
        self.pref_rad_py = np.array(self.pref_rad, dtype=np.int32, order='F')
        self.numlp_py = np.array(self.numlp, dtype=np.int32, order='F')
        self.numradlp_py = np.array(self.numradlp, dtype=np.int32, order='F')

        print("nstrt1",self.nstrt1,self.strt_struc1)

    def get_orbs(self, orb_class):
        self.atoset, self.norbsym, self.active, self.atn, self.orbsym, self.actv_atm_num\
                = orb_class.get_orbital_matrices()
        print('orbs_datai,atoset', self.atoset)
        print('norbsym', self.norbsym)
        print('active', self.active)
        print('atn', self.atn)
        print('orbsym', self.orbsym)
        print('actv_atm_num', self.actv_atm_num)

    def get_ctrl_keywds(self, ctrl_keywds):
        self.geometry_unit, self.nao, self.nae, self.nmul = ctrl_keywds.get_ctrl_keywds()
        self.nao_py = int(self.nao)
        self.nae_py = int(self.nae)
        self.nmul_py = int(self.nmul)

    def get_geometry(self):
        if not GlobVar.geometry_inserted:
            messagebox.showerror("Invalid Geometry", "please insert the geometry")
            return

        symat, coordx, coordy, coordz, symatno = GlobVar.readgeo.get_geometry_data()

        #bohrtoangs = 1.0
        #if GlobVar.geo_unit == 'Bohr':
        #    bohrtoangs = 0.529177

        # Initialize fixed-size arrays
        self.symat_py = np.zeros(GlobVar.total_atoms, dtype="U5", order='F')
        self.coordx_py = np.zeros(GlobVar.total_atoms, dtype=np.float64, order='F')
        self.coordy_py = np.zeros(GlobVar.total_atoms, dtype=np.float64, order='F')
        self.coordz_py = np.zeros(GlobVar.total_atoms, dtype=np.float64, order='F')
        self.symatno_py = np.zeros(GlobVar.total_atoms, dtype=np.float64, order='F')

        #self.nstrt1 = 2

        print("strt_struc1",self.strt_struc1)
        print("main_bond_py",self.main_bond_py)
        print("pref_rad_lp_py",self.pref_rad_lp_py)
        print("pref_rad_py",self.pref_rad_py)
        print("numlp_py",self.numlp_py)
        print("numradlp_py",self.numradlp_py)
        print("numrad_py",self.numrad_py)
        print("numudr_py",self.numudr_py)
        #self.symat_py = np.zeros(100, dtype="S5", order='F')
        #self.coordx_py = np.zeros(100, dtype=np.float64, order='F')
        #self.coordy_py = np.zeros(100, dtype=np.float64, order='F')
        #self.coordz_py = np.zeros(100, dtype=np.float64, order='F')
        #self.symatno_py = np.zeros(100, dtype=np.float64, order='F')

        # Number of atoms
        #n = len(coordx)  

        # Copy and convert data with scaling
        #self.symat_py[:n] = symat
        #self.coordx_py[:n] = np.array(coordx) #* bohrtoangs
        #self.coordy_py[:n] = np.array(coordy) #* bohrtoangs
        #self.coordz_py[:n] = np.array(coordz) #* bohrtoangs
        #self.symatno_py[:n] = np.array(symatno, dtype=np.float64)

        self.symat_py[:] = np.asarray(symat, dtype='U5')
        self.coordx_py[:] = np.asarray(coordx, dtype=np.float64)
        self.coordy_py[:] = np.asarray(coordy, dtype=np.float64)
        self.coordz_py[:] = np.asarray(coordz, dtype=np.float64)
        self.symatno_py[:] = np.asarray(symatno, dtype=np.float64)

#       # print('geometry:',self.symat_py, self.coordx_py, self.coordy_py, self.coordz_py, self.symatno_py)


    def share_input_data(self, output_folder):
        self.out_folder_name = output_folder
        print(self.norbsym)
        print(self.geometry_unit)
        print(self.nao_py)
        print(self.nae_py)
        print(self.nmul_py)
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
        print(GlobVar.total_atoms)
        print(GlobVar.num_iao)
        print(self.actv_atm_num)
        print(self.out_folder_name)
        print(self.strt_struc)
        print(self.main_bond_py)
        print(self.pref_rad_lp_py)
        print(self.pref_rad_py)
        print(self.numlp_py)
        print(self.numrad_py)
        print(self.numradlp_py)
        print(self.numudr_py)
        symm_str.ctrl_mod.get_ctrl_inputs(
                self.nao_py,
                self.nae_py,
                self.nstrt1,
                GlobVar.total_atoms,
                self.geometry_unit,
                self.nmul_py,
                self.chinst, 
                self.symm,
                self.asymm,
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
                GlobVar.num_iao,
                self.actv_atm_num,
                self.strt_struc,
                self.main_bond_py,
                self.pref_rad_lp_py,
                self.pref_rad_py,
                self.numlp_py,
                self.numrad_py,
                self.numradlp_py,
                self.numudr_py,
                GlobVar.molecule_string,
                self.out_folder_name
                )
        return (self.mout)


class Output:
    def __init__(self, root):
        self.root = root

    def load_structure_file(self, fname):
        self.fname = fname
        Output_window = tk.Toplevel(self.root)
        Output_window.title("Output")
        Output_window.geometry("1050x950")
        Output_window.configure(background="lightblue")
        structures_cov = []
        various_qualities_cov = []
        overall_qualities_cov = []
        rumers_cov = []
        sls_cov = []
        structures_ion = []
        various_qualities_ion = []
        overall_qualities_ion = []
        rumers_ion = []
        sls_ion = []
        sets= perm_cov_str= perm_ion_str= all_cov_str= all_ion_str = 0
        print('I am here')

        self.frames = {}

        frame_info = {
                "frame":("lightblue", 0, 0),
                "text_frame":("white",1, 0),
                "button_frame":("lightblue",2,0),
                "set_frame":("white",3,0),
                "info_frame":("white",3,1),
                "info_button_frame":("lightblue",4,0),
                }
        for frame_name,(bg, row, column) in frame_info.items():
            frame = tk.Frame(Output_window, bg=bg)
            frame.grid(row=row, column=column, padx=10, pady=10)
            self.frames[frame_name] = frame

        self.frame_view_str = tk.Frame(self.frames['info_button_frame'], bg='lightblue')
        self.frame_view_str.grid(row=0, column=0,padx=20, pady=20)
        self.frame_view_info = tk.Frame(self.frames['info_button_frame'], bg='lightblue')
        self.frame_view_info.grid(row=0, column=1,padx=20, pady=20)
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
        # calculation number of available structures & allowed structures
        nao = GlobVar.num_orbital
        nae = GlobVar.num_electron
        nlp =abs(nao - nae)
        nmult = GlobVar.multiplicity
        cov_nlp = nlp
        nlast = nmult - 1
        cov_bond = int((nao - 2 * cov_nlp - nlast)/2)
        if GlobVar.CovIon==1 or GlobVar.CovIon==2:
            sets, perm_cov_str, all_cov_str = self.wigner(nlp, nao, nae, nmult)
        if GlobVar.CovIon==1 or GlobVar.CovIon==3:
            for i in range(1, cov_bond+1):
                ion_nlp = cov_nlp + i
                w1 , w2, w3 = self.wigner(ion_nlp, nao, nae, nmult)
                perm_ion_str += w2
                all_ion_str += w3 

        filename = fname + '/' + 'structures.dat'

        if not os.path.exists(filename):
            text_widget.insert(tk.END, "No Structures are available\n")
            return
        else:
            sw = GlobVar.num_orbital*2+20
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
            i = 0
            with open(filename, "r") as file:
                lines = file.readlines()
                for line in lines:
                    if line.startswith('Cov Space'):
                        cols = line.strip().split()
                        space_num = cols[3]
                        result = (f"Cov Space = {space_num}\n\n")
                        covionflg = 1
                        sl = 0
                        structures = []
                        various_qualities = []
                        overall_qualities = []
                        rumers = []
                        sls = []
                        text_widget.insert(tk.END, result)
                    if line.startswith('Ion Space'):
                        cols = line.strip().split()
                        space_num = cols[3]
                        result = (f"Ion Space = {space_num}\n\n")
                        covionflg = 2
                        #sl = 0
                        structures = []
                        various_qualities = []
                        overall_qualities = []
                        rumers = []
                        sls = []
                        text_widget.insert(tk.END, result)
                    if line.startswith("structure"):
                        # print('line', line)
                        cols = line.strip().split()
                        if covionflg == 1:
                            sl = cols[1]
                            sls.append(int(sl))
                            various_quality = " ".join(cols[3:8])
                            various_qualities.append(various_quality)
                            print('various_quality',various_quality)
                            overall_quality = cols[10]
                            overall_qualities.append(int(overall_quality))
                            rumer = cols[12]
                            rumers.append(rumer)
                            if rumer == 'Rumer':
                                rumer = 'R'
                            structure = " ".join(cols[13:])
                            structures.append(structure)
                            sw = len(structure) + 10
                            #i += 1
                            result = (
                                    f"{sl:<16}"
                                    f"{structure:^{sw}}"
                                    f"{various_quality:^30}"
                                    f"{overall_quality:^22}"
                                    f"{rumer:^20}\n\n"
                                    )
                            text_widget.insert(tk.END, result)
                        if covionflg == 2:
                            sl = cols[1]
                            sls.append(int(sl))
                            various_quality = " ".join(cols[3:8])
                            various_qualities.append(various_quality)
                            print('various_quality',various_quality)
                            overall_quality = cols[10]
                            overall_qualities.append(int(overall_quality))
                            rumer = cols[12]
                            rumers.append(rumer)
                            if rumer == 'Rumer':
                                rumer = 'R'
                            structure = " ".join(cols[13:])
                            structures.append(structure)
                            sw = len(structure) + 10
                            #i += 1
                            result = (
                                    f"{sl:<16}"
                                    f"{structure:^{sw}}"
                                    f"{various_quality:^30}"
                                    f"{overall_quality:^22}"
                                    f"{rumer:^20}\n\n"
                                    )
                            text_widget.insert(tk.END, result)
                    if line.startswith(" ============================================================"):
                        if covionflg == 1:
                            structures_cov.append(structures)
                            various_qualities_cov.append(various_qualities)
                            overall_qualities_cov.append(overall_qualities)
                            rumers_cov.append(rumers)
                            sls_cov.append(sls)
                        if covionflg == 2:
                            structures_ion.append(structures)
                            various_qualities_ion.append(various_qualities)
                            overall_qualities_ion.append(overall_qualities)
                            rumers_ion.append(rumers)
                            sls_ion.append(sls)
                print('structures_cov',structures_cov)
                print('structures_ion',structures_ion)

                label_tot_cov_str = ttk.Label(
                        self.frames["frame"], text=(f'Number of available covalent structures of thes ystem = {all_cov_str}\n'
                                                    f"Number of available ionic structures of the system = {all_ion_str}"), 
                        style="Colour_Label.TLabel"
                        )
                label_tot_cov_str.grid(row=0, column=0, padx=5, pady=5)
                label_tot_all_cov_str = ttk.Label(
                        self.frames["frame"], text=(f'Number of allowed covalent structures = {perm_cov_str}\n'
                                                    f'Number of allowed ionic structures = {perm_ion_str}'), 
                        style="Colour_Label.TLabel"
                        )
                label_tot_all_cov_str.grid(row=0, column=1, padx=5, pady=5)
        text_widget.config(state=tk.DISABLED)

        self.tempfname = fname + '/' + 'out.temp'

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

        print('out.temp file path', self.tempfname)
        if not os.path.exists(self.tempfname):
            messagebox.showerror(
                    'Error', 'No Output File has been found'
                    )
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
###################################################################################################

       ## Below part is to find out highest quality structures and update the list(lines) ##
       ## for the option 'Best_Sets'.

        if GlobVar.num_sets == 1:
            if GlobVar.CovIon == 1 or GlobVar.CovIon == 2:
                best_strc_cov = []
                strc_cov = None

            if GlobVar.CovIon == 1 or GlobVar.CovIon == 3:
                best_strc_ion = []
                strc_ion = None

            for line in lines:
                if GlobVar.CovIon == 1 or GlobVar.CovIon == 2:
                    if line.startswith('Cov_Space'):
                        if strc_cov is not None:
                            best_strc_cov.append(strc_cov)
                        strc_cov = []
                    if line.startswith('Set_number_c'):
                        linear_indset = line.strip().split()
                        strc_cov.append(linear_indset[2:])

                if GlobVar.CovIon == 1 or GlobVar.CovIon == 3:
                    if line.startswith('Ion_Space'):
                        if strc_ion is not None:
                            best_strc_ion.append(strc_ion)
                        strc_ion = []
                    if line.startswith('Set_number_i'):
                        linear_indset = line.strip().split()
                        strc_ion.append(linear_indset[2:])

            if GlobVar.CovIon == 1 or GlobVar.CovIon == 2:
                if strc_cov is not None:
                    best_strc_cov.append(strc_cov)
         #       print('best_strc_cov',best_strc_cov)

            if GlobVar.CovIon == 1 or GlobVar.CovIon == 3:
                if strc_ion is not None:
                    best_strc_ion.append(strc_ion)
          #      print('best_strc_ion',best_strc_ion)
          #  print('overall_qualities_cov',overall_qualities_cov)
          #  print('overall_qualities_ion',overall_qualities_ion)

            

            if GlobVar.CovIon == 1 or GlobVar.CovIon == 2:
                result_cov = []
                for struct, qualities in zip(best_strc_cov, overall_qualities_cov):
                
                    mapped_sums = []
                    min_indices = []
                
                    for group in struct:  
                        values = [qualities[int(x)-1] for x in group]
                        mapped_sums.append(sum(values))
                
                    if mapped_sums:
                        min_val = min(mapped_sums)
                        min_indices = [i for i, x in enumerate(mapped_sums) if x == min_val]
                        min_indices = [i + 1 for i in min_indices]
                    else:
                        normalized = []
                
                    #result_cov.append(normalized)
                    result_cov.append(min_indices)
            
            if GlobVar.CovIon == 1 or GlobVar.CovIon == 3:
                result_ion = []
                for struct, qualities in zip(best_strc_ion, overall_qualities_ion):
                
                    mapped_sums = []
                    min_indices = []
                
                    for group in struct:  
                        values = [qualities[int(x)-1] for x in group]
                        mapped_sums.append(sum(values))
                
                    if mapped_sums:
                        min_val = min(mapped_sums)
                        min_indices = [i for i, x in enumerate(mapped_sums) if x == min_val]
                        min_indices = [i + 1 for i in min_indices]
                    else:
                        normalized = []
                
                    result_ion.append(min_indices)
            
        #    print('result_cov',result_cov)
        #    print('result_cov_1',result_cov_1)
            cov_block_idx = 0
            ion_block_idx = 0
            
            lines1 = [] # temporary list to store lowest structure indices
            i = 0
            while i < len(lines):
            
                line = lines[i]
                #print('line',line)
            
                # -------- Cov_Space --------
                if (GlobVar.CovIon == 1 or GlobVar.CovIon == 2) and line.startswith('Cov_Space'):
                    #if line.startswith('Cov_Space'):
                    lines1.append(line)  
                    #print('lines1',lines1)
                        
                    i += 1
                    set_idx = 1  
                    counter = 1
            
                    while i < len(lines) and lines[i].startswith('Set_number_c'):
                        if cov_block_idx < len(result_cov) and set_idx in result_cov[cov_block_idx]:
                            parts = lines[i].split()
                            #print('parts',parts)

                            new_line = f"Set_number_c  {counter}   " + " ".join(parts[2:]) + "\n"
                            #print('new_line',new_line)
                            lines1.append(new_line) 
                            #print('lines1',lines1)    
                            counter += 1
                            
                        set_idx += 1
                        i += 1
            
                    cov_block_idx += 1
                    continue

                # -------- Ion_Space --------
                if (GlobVar.CovIon == 1 or GlobVar.CovIon == 3) and line.startswith('Ion_Space'):
                    #if line.startswith('Ion_Space'):
                    lines1.append(line)  # keep header

                    i += 1
                    set_idx = 1
                    counter = 1

                    while i < len(lines) and lines[i].startswith('Set_number_i'):

                        if ion_block_idx < len(result_ion) and set_idx in result_ion[ion_block_idx]:
                            parts = lines[i].split()

                            new_line = f"Set_number_i  {counter}   " + " ".join(parts[2:]) + "\n"
                            lines1.append(new_line)

                            #print('lines1*',lines1)    
                            counter += 1

                        set_idx += 1
                        i += 1

                    ion_block_idx += 1
                    continue
                #sys.exit()
            lines = lines1.copy()
        #print('lines1****',lines)
###########################################################################################

        flag1 = 1
        flag2 = 2
        flag3 = 3
        view_set_button_cov = tk.Button(
                self.frames["button_frame"], 
                text='View Cov Set', 
                command=lambda: self.view_set(flag1, structures_cov, various_qualities_cov, overall_qualities_cov,
                                              rumers_cov, sls_cov, set_text_wt, lines, 'cov')
                )
        view_set_button_cov.grid(row=0, column=2, padx=10, pady=10)
        if GlobVar.CovIon == 3:
            view_set_button_cov.config(state=tk.DISABLED)  # Initially disable the button

        view_set_button_ion = tk.Button(
                self.frames["button_frame"], 
                text='View Ion Set', 
                command=lambda: self.view_set(flag1, structures_ion, various_qualities_ion, overall_qualities_ion,
                                              rumers_ion, sls_ion, set_text_wt, lines, 'ion')
                )
        view_set_button_ion.grid(row=0, column=3, padx=10, pady=10)
        if GlobVar.CovIon == 2:
            view_set_button_ion.config(state=tk.DISABLED)  # Initially disable the button
        
        def handle_next_set():
            if GlobVar.ciflg == 'cov':
                self.view_set(flag2, structures_cov, various_qualities_cov, overall_qualities_cov,
                            rumers_cov, sls_cov, set_text_wt, lines, '')
            elif GlobVar.ciflg == 'ion':
                self.view_set(flag2, structures_ion, various_qualities_ion, overall_qualities_ion,
                            rumers_ion, sls_ion, set_text_wt, lines, '')

        def handle_prev_set():
            if GlobVar.ciflg == 'cov':
                self.view_set(flag3, structures_cov, various_qualities_cov, overall_qualities_cov,
                            rumers_cov, sls_cov, set_text_wt, lines, '')
            elif GlobVar.ciflg == 'ion':
                self.view_set(flag3, structures_ion, various_qualities_ion, overall_qualities_ion,
                            rumers_ion, sls_ion, set_text_wt, lines, '')


    # Create the button and assign the nested function
        view_prevtset_button = tk.Button(
            self.frames["button_frame"],
            text='Prev Set',
            command=handle_prev_set
        )
        view_prevtset_button.grid(row=0, column=1, padx=10, pady=10)
    #CovIon = None
        if GlobVar.num_sets == 0:
            view_prevtset_button.config(state=tk.DISABLED)  # Initially disable the button


        view_nextset_button = tk.Button(
            self.frames["button_frame"], 
            text='Next Set', 
            command= handle_next_set
        )
        view_nextset_button.grid(row=0, column=4, padx=10, pady=10)

        if GlobVar.num_sets == 0:
            view_nextset_button.config(state=tk.DISABLED)  # Initially disable the button


    def view_set(self, flag, structure_ci, various_qualities_ci, overall_qualities_ci, 
                 rumers_ci, sls_ci, set_text_wt, lines, covionflg):

        various_qualities_selected = []
        overall_qualities_selected = []
        rumers_selected = []
        indset = []
        diff_sets = 0
        max_set = 0

        if covionflg != '':
            GlobVar.ciflg = covionflg

        #print('ciflg********',GlobVar.ciflg)

        if flag == 1:
            GlobVar.set_id = 1
            self.set_info_wt.delete("1.0", "end")
        elif flag == 2:
            GlobVar.set_id += 1
            self.set_info_wt.delete("1.0", "end")
        elif flag == 3:
            GlobVar.set_id -= 1
            self.set_info_wt.delete("1.0", "end")

        set_text_wt.tag_configure("right", justify="right")
        count_c = 0
        count_i = 0

        for line in lines:
            if line.startswith('Set_number_c') and GlobVar.ciflg == 'cov':
                linear_indset = line.strip().split()
                diff_sets = linear_indset[1]
            if line.startswith('Set_number_i') and GlobVar.ciflg == 'ion':
                linear_indset = line.strip().split()
                diff_sets = linear_indset[2]

        max_set = diff_sets
        linear_indset = 0
        structure = []
        various_qualities=[]
        overall_qualities=[]
        rumers=[]
        sls=[]

        for line in lines:
            if 'Cov_Space' in line and GlobVar.ciflg == 'cov':
                structure=structure_ci[count_c]
                various_qualities=various_qualities_ci[count_c] 
                overall_qualities=overall_qualities_ci[count_c] 
                rumers=rumers_ci[count_c]
                sls=sls_ci[count_c]
                count_c += 1
            if 'Ion_Space' in line and GlobVar.ciflg == 'ion':
                structure=structure_ci[count_i]
                various_qualities=various_qualities_ci[count_i] 
                overall_qualities=overall_qualities_ci[count_i] 
                rumers=rumers_ci[count_i]
                sls=sls_ci[count_i]
                count_i += 1
            if line.startswith('Set_number_c') and GlobVar.ciflg == 'cov':
                linear_indset = line.strip().split()

                if linear_indset[1] == str(GlobVar.set_id):
                    if count_c == 1:
                        #indset = []
                        set_text_wt.delete("1.0", "end")
                        set_text_wt.insert(
                                tk.END, f"\n    Total Number of  Independent Set of Structures: {max_set}  \n\n"
                                )
                        set_text_wt.insert(
                                tk.END, "\n           Set number::  " + f"{linear_indset[1]}\n\n"
                                )
                    sl = 0
                    for index in linear_indset[2:]:
                        #print('linear_indset',linear_indset[2:])
                        #print('structure',structure[int(index) - 1],[int(index) - 1])
                        sl += 1
                        indset.append(int(index))
                        set_text_wt.insert(
                                tk.END, f" {sl}  ({sls[int(index) - 1]}):     " + structure[int(index) - 1] + "\n"
                                )
                        various_qualities_selected.append(various_qualities[int(index) - 1])
                        overall_qualities_selected.append(overall_qualities[int(index) - 1])
                        self.set_info(various_qualities_selected, overall_qualities_selected)

                    #info_button = tk.Button(self.frame_view_info, text="Set Info", 
                    #       command=lambda: self.set_info(
                    #            various_qualities_selected, overall_qualities_selected
                    #            )
                    #        )
                    #info_button.grid(row=0, column=0, sticky=tk.W, padx=10)
                    if GlobVar.symm_key == 1:
                        symm_grp_button = tk.Button(self.frame_view_info, text="Symmetry groups",
                                                    command=self.show_symm_groups
                                                    )
                        symm_grp_button.grid(row=0, column=1, sticky=tk.W, padx=10)

                    ## Structure view Buttons ... ##
                    self.molkey = 0

                    view_button = tk.Button(self.frame_view_str, text="View Structure", 
                                            command=lambda i=indset, s=structure, l=linear_indset[1]:
                                            self.View_Structure(i, 'molflg', s, l))
                    view_button.grid(row=0, column=1, sticky=tk.W, padx=10)

                    view_button_next = tk.Button(self.frame_view_str, text="next", 
                                            command=lambda i=indset, s=structure, l=linear_indset[1]:
                                            self.View_Structure(i, 'molflg1', s, l))
                    view_button_next.grid(row=0, column=2, sticky=tk.W, padx=10)

                    view_button_prev = tk.Button(self.frame_view_str, text="prev", 
                                            command=lambda i=indset, s=structure, l=linear_indset[1]:
                                            self.View_Structure(i, 'molflg2', s, l))
                    view_button_prev.grid(row=0, column=0, sticky=tk.W, padx=10)
                    #print('linear_indset', linear_indset)

            if line.startswith('Set_number_i') and GlobVar.ciflg == 'ion':
                linear_indset = line.strip().split()
                print('linear_indset',linear_indset)

                if linear_indset[1] == str(GlobVar.set_id):
                    if count_i == 1:
                        set_text_wt.delete("1.0", "end")
                        set_text_wt.insert(
                                tk.END, f"\n    Total Number of  Independent Set of Structures: {max_set}  \n\n"
                                )
                        set_text_wt.insert(
                                tk.END, "\n           Set number::  " + f"{linear_indset[1]}\n\n"
                                )
                    sl = 0
                    for index in linear_indset[2:]:
                        sl += 1
                        indset.append(int(index))
                        set_text_wt.insert(
                                tk.END, f" {sl}  ({sls[int(index) - 1]}):        " + structure[int(index) - 1] + "\n"
                                )
                        various_qualities_selected.append(various_qualities[int(index) - 1])
                        overall_qualities_selected.append(overall_qualities[int(index) - 1])
                        self.set_info(various_qualities_selected, overall_qualities_selected)

                    #info_button = tk.Button(self.frame_view_info, text="Set Info", 
                    #       command=lambda: self.set_info(
                    #            various_qualities_selected, overall_qualities_selected
                    #            )
                    #        )
                    #info_button.grid(row=0, column=0, sticky=tk.W, padx=10)
                    if GlobVar.symm_key == 1:
                        symm_grp_button = tk.Button(self.frame_view_info, text="Symmetry groups",
                                                    command=self.show_symm_groups
                                                    )
                        symm_grp_button.grid(row=0, column=1, sticky=tk.W, padx=10)

                    ## Structure view Buttons ... ##
                    self.molkey = 0

                    view_button = tk.Button(self.frame_view_str, text="View Structure", 
                                            command=lambda i=indset, s=structure, l=linear_indset[1]:
                                            self.View_Structure(i, 'molflg', s, l))
                    view_button.grid(row=0, column=1, sticky=tk.W, padx=10)

                    view_button_next = tk.Button(self.frame_view_str, text="next", 
                                            command=lambda i=indset, s=structure, l=linear_indset[1]:
                                            self.View_Structure(i, 'molflg1', s, l))
                    view_button_next.grid(row=0, column=2, sticky=tk.W, padx=10)

                    view_button_prev = tk.Button(self.frame_view_str, text="prev", 
                                            command=lambda i=indset, s=structure, l=linear_indset[1]:
                                            self.View_Structure(i, 'molflg2', s, l))
                    view_button_prev.grid(row=0, column=0, sticky=tk.W, padx=10)
                    #print('linear_indset', linear_indset)

    def View_Structure(self, indset, flg, structures, set_number):
        print('indset',indset, flg)
        print('structure',structures[0], structures[1])
        seen = set()
        nlp_lst = []
        rad_lst = []
        nbd = []
        bond_lst = []

        if flg == 'molflg':
            self.molkey = 0
        if flg == 'molflg1':
            self.molkey += 1
        if flg == 'molflg2' and self.molkey != 0:
            self.molkey -= 1
        if self.molkey < len(indset):
            entry = indset[self.molkey] - 1
            structure =structures[entry].split() 
            print('entry', entry)
            print(structures[entry])

            nlp = abs(GlobVar.num_electron - GlobVar.num_orbital)
            if nlp != 0:
                for item in structure:
                    if item not in seen:
                        seen.add(item)
                        nlp_lst.append(item)
                    if len(result) == nlp:
                        break

            nrad = GlobVar.multiplicity - 1
            if nrad != 0:
                rad_lst = structure[-nrad:]

            nbond = int((len(structure) - nrad - 2*nlp)/2)
            if nbond != 0:
                k=nrad
                usable = structure[:-k] if k > 0 else structure
                bond_lst = [usable[i:i+2] for i in range(len(usable) - 2*nbond, len(usable), 2)]
        else:
            return
        d=drawing_molecule()
        d.str_info(nlp_lst, rad_lst, bond_lst, self.molkey+1, set_number)
        d.VB_view()

        print('nlp:',nlp_lst,'rad:',rad_lst,'bond:',bond_lst)


    def show_symm_groups(self):
        top = tk.Toplevel(root)   # create new window
        top.title("Symmetry sub-groups")
        top.geometry("500x600")
        top.configure(background="lightblue")

        canvas = tk.Canvas(top, bg="lightblue")
        scrollbar = tk.Scrollbar(top, orient="vertical", command=canvas.yview)
        scroll_frame = tk.Frame(canvas)
        scroll_frame.configure(background="lightblue")

        scroll_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
            )

        #canvas = tk.Canvas(container, height=370, width=510)
        canvas.create_window((0, 0), window=scroll_frame,  anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")


        input_file = self.fname + '/' + 'out_symm.temp'
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
            current_space = None
            current_block = []
            results = []
        
            for line in lines:
                line = line.strip()
        
                # Detect new space
                if line.startswith("Cov Space") or line.startswith("Ion Space"):
                    if current_space is not None:
                        results.append((current_space, current_block))
                    current_space = line
                    current_block = []
                    continue
        
                # Skip empty lines
                if not line:
                    continue
        
                parts = line.split()
        
                # Skip malformed lines
                if len(parts) < 3:
                    continue
        
                serial = int(parts[0])
                quality = int(parts[1])
                structure = parts[2:]
        
                current_block.append((serial, quality, structure))
        
            # Append last block
            if current_space is not None:
                results.append((current_space, current_block))
        grouped_results = []

        for space, block in results:
            groups = defaultdict(list)

            for serial, quality, structure in block:
                groups[quality].append((serial, structure))

            grouped_results.append((space, dict(groups)))
        
        for space, groups in grouped_results:

            tk.Label(scroll_frame,
                     text=space,
                     font=("Arial", 16, "bold"),
                     bg="lightblue",
                     anchor="w").pack(fill="x", padx=10, pady=5)

            for q in sorted(groups.keys()):

                tk.Label(scroll_frame,
                         text=f"  Sub-group (Symmetry score {q}):",
                         font=("Arial", 14, "bold"),
                         bg="lightblue",
                         anchor="w").pack(fill="x", padx=20)

                for serial, structure in groups[q]:
                    struct_str = " ".join(structure)

                    tk.Label(scroll_frame,
                            text=f"    {serial}:  {struct_str}",
                            font=("Arial", 14),
                            bg="lightblue",
                            anchor="w",
                            justify="left").pack(fill="x", padx=30)

            tk.Label(scroll_frame, text="", bg="lightblue").pack(pady=5)

        tk.Button(top, text="Close", command=top.destroy).pack(pady=5)

    def set_info(self, various_qualities, overall_qualities):
        oqt = sum(overall_qualities[i] for i in range(len(overall_qualities)))
        self.set_info_wt.delete("1.0", "end")
        self.set_info_wt.insert(tk.END, f"\n Overall Quality: {oqt}\n\n")
        self.set_info_wt.insert(tk.END, " IAB  NAB  SBB  PDB  PDR\n")
        self.set_info_wt.insert(tk.END, f"{' -----':<8}" 
                                f"{'-----':<8}" 
                                f"{'-----':<8}" 
                                f"{'-----':<8}" 
                                f"{'-----':<8}\n")
        for i in range(len(various_qualities)):
            a, b, c, d, e = various_qualities[i].strip().split()
            self.set_info_wt.insert(tk.END,f"  {a:<8}" 
                                    f"{b:<8}" 
                                    f"{c:<8}" 
                                    f"{d:<8}" 
                                    f"{e:<8}""\n")
        if GlobVar.Rum_Ch == 0:
            self.set_info_wt.insert(tk.END, "\n The set is a Rumer Set\n")
        elif GlobVar.Rum_Ch == 1:
            self.set_info_wt.insert(tk.END, "\n The set is a Chemical Insight Set \n")
        elif GlobVar.Rum_Ch == 2:
            self.set_info_wt.insert(tk.END, "\n The set is a Equal Bond-Distributed Set \n")

        return #oqt

    def wigner(self, nlp, nao, nae, nmult):
        ''' none = number of one electron orbital
            nlp = number of lone paires
            nao = number of active electrons
            nmult = multiplicity
        '''
        print("nlpppp",nlp)
        wigner = 0
        sets = 1

        ''' calculate permissible number of structre with Wigner's theorem'''
        noeo = nae - 2*nlp
        navo = nao - nlp
        c = math.comb(nao, navo)
        c1 = math.comb(navo, noeo)
        spin = (nmult - 1)/2
        neum = nmult*math.factorial(noeo)
        denom1 = math.factorial(int((noeo/2) + spin+1))
        denom2 = math.factorial(int((noeo/2) - spin))
        wigner = int(neum/(denom1 * denom2))*c*c1


        ''' Number of sets depending on lone pairs'''
        print("nlp*****",nlp)
        if nlp != 0:
            sets = math.comb(nao, nlp)

        ''' Calculate total number of structures'''
        term1 = math.comb(nao, nlp)
        term2 = math.comb(noeo, int(2*spin))
        term3 = math.comb(navo, noeo)

        product_term = np.prod([(math.comb(noeo - int(2*spin) - 2*i, 2))
                                for i in range(int((noeo / 2)-spin))])
        denom = math.factorial(int((noeo / 2) - spin))

        totstr = int((term1 * term2 * term3 * product_term) / denom)

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
    root.geometry("630x950")
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
