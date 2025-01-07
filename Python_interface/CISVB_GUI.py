################################################################################################################
# python GUI interface to interact with the users to get input keywords geometry file orbital informations and #
# send them to Fortran subroutines to generate set of VB structures and store them in the output file.         #
################################################################################################################
    
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
import cisvb
import numpy as np
import Pmw

num_orbital, num_electron, multiplicity = None, None, None
type_orb_count = None
num_iao = None
geometry_inserted = True
at_list_bold = [
    'H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE', 'NA', 'MG', 'AL',
    'SL', 'P', 'S', 'CL', 'AR', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NR',
    'CU', 'ZN', 'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC',
    'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', 'CE',
    'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA',
    'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA'
]


class Ctrl_Input:
    def __init__(self, root):
        self.ctrl_inputs={}   # initialise disctionary to store control data
        self.root = root
        self.input_text=''
        self.insert = False
        self.orbital_button = None
        self.structure_type_entry = False
        self.method_type_entry = False
        self.unit_type_entry = False
        self.readgeo = None
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
        lebel = ttk.Label(self.frame1, text = "Brows to Upload Geometry File", style = "Colour_Label.TLabel")
        lebel.grid(row = 0, column=0, sticky=tk.W, padx=5, pady=5)
        button = ttk.Button(self.frame1, text = "Brows", command = self.read_geometry)
        button.grid(row = 0, column = 1, sticky = tk.W, padx = 5, pady = 5)

        # upload geometry manually
        lebel = ttk.Label(self.frame1, text = "Insert Geometry Manually", style = "Colour_Label.TLabel")
        lebel.grid(row = 2, column=0, sticky=tk.W, padx=5, pady=5)
        button = ttk.Button(self.frame1, text = "Geometry", command = self.insert_geo_manually)
        button.grid(row = 2, column = 1, sticky = tk.W, padx = 5, pady = 5)

        insert_button = ttk.Button(self.frame4, text="Insert", command=self.validate_and_generate)
        insert_button.grid(row=4, column=1, padx=5, pady=10)
        self.balloon.bind(insert_button, "Press Insert after providing all control keywords")

        self.tip_label = ttk.Label(self.frame3, text="Select structure type: ", style = "Colour_Label.TLabel")
        self.tip_label.grid(row = 0, column = 0, pady=10)
        self.balloon.bind(self.tip_label,"Select any one from the three types. select 'Covalent'\n"
                                         " to generate covalent sets only. Select 'Ionic' to generate \n" 
                                         "ionic structures sets only. Select both to generate both \n"
                                          "type of structure sets.")
        self.tip_method_label = ttk.Label(self.frame3, text="Select method type: ", style = "Colour_Label.TLabel")
        self.tip_method_label.grid(row = 1, column = 0, pady=10)
        self.balloon.bind(self.tip_method_label,"Select any one from the two method provided. select 'Rumer' to genarate a Rumer structure\
set, no need to provide any spatial keywords. For all other type of set please select Chem Inst.")

        # Variable to store the selected input
        self.str_type = tk.StringVar(value="None")
        self.method_type = tk.StringVar(value="None")

        # Create radio buttons
        self.str_tip_buttons()

    def create_geo_unit(self):
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

#        if self.file_name == None:
#            self.file_name='CISVB_OUT.DAT'
#        return (self.file_name)

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
                if element in at_list_bold:
                    atomic_number = at_list_bold.index(element) + 1  # Atomic number = index + 1
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
        self.file_path = filedialog.askopenfilename(title="Select Geometry File")

        if self.file_path:  # Check if a file is selected
            readgeo = Read_Geo(self.file_path)  # Initialize the Read_Geo class
            readgeo.read_geometry()       # Call the read_geometry method

            ttk.Button(
                self.frame1,
                text=f"View_Geometry",
                command=self.display_geometry_file
            ).grid(row=0, column=2, columnspan=2)



    def display_geometry_file(self):
        with open(self.file_path, "r") as file:
            content = file.read()

        # Create a Toplevel window to display the file content
        geo_file_display = tk.Toplevel()
        geo_file_display.title("geometry File Content")
        geo_file_display.geometry("600x400")  # Optional: Set the size of the Toplevel window

        # Add a Text widget to display file content
        text_widget = tk.Text(geo_file_display, wrap="word")
        text_widget.insert("1.0", content)  # Insert content into the Text widget
        text_widget.pack(expand=True, fill="both")



    def insert_geo_manually(self):
        """Allows manual insertion of geometry data via Read_Geo."""
        self.readgeo = Read_Geo(self.file_path)
        self.readgeo.insert_geo(self.root)  # Call insert_geo method
        if geometry_inserted == True:
            button = ttk.Button(self.frame1, text = "View Geometry", command = self.readgeo.display_geometry)
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
        global num_iao
        self.ctrl_inputs = []
        try:
            molecule_string = self.molecule_entry.get()
            Total_Electrons = self.count_Total_Electron(molecule_string)
            print('Total Number of Electrons', Total_Electrons)
        except ValueError:
            massagebox.showerror("Molecule Error","Please enter molecular Formula Properly ")


        if self.structure_type_entry == False:
            messagebox.showerror("Structure Type Missing", "Please Select a Structure Type among Covalent, Ionic, and Both.")
            return
        if self.method_type_entry == False:
            messagebox.showerror("Method Type Missing", "Please Select a Method between Chemical Insight and Rumer.")
            return
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
                print('nao, nae, nmul',nao, nae, nmul)
                if (Total_Electrons- int(nae)) % 2== 0:
                    num_iao = int((Total_Electrons- int(nae))/2)
                    print('num_iao',num_iao)
                else:
                    massagebox.showerror("Active Orbital Error","The number of Active orbitals or \n"
                                         "Molecular Formula is not correct")
                self.orbital_button.config(state=tk.NORMAL)
            return True  # Validation successful

    def generate_ctrl_input(self):
        if not self.insert:
            messagebox.showerror("Validation Error", "Please insert valid values before insert other data.")
            return
        nao, nae, nmul = self.ctrl_inputs
        # Collect all inputs
#        print("Control Inputs:", self.ctrl_inputs)
#        for key, entry in self.entries.items():
#            self.ctrl_inputs[key] = entry.get().strip()
        return(nao, nae, nmul)

    def str_tip_buttons(self):
        # Tips to display
        tips = ["Covalent", "Ionic", "Both"]
        tip_method = ["Chem inst", "Rumer"]

        # Create and pack each radio button
        i=2
        for tip in tips:
            i=i+1
            button = ttk.Radiobutton(
                self.frame3,
                text=tip,
                value=tip,
                variable=self.str_type,
                command=self.update_str_type,
                style="Custom.TRadiobutton"
            )
            button.grid(row = 0, column = i, padx=10, pady=10)

        i=2
        for tip in tip_method:
            i=i+1
            button = ttk.Radiobutton(
                self.frame3,
                text=tip,
                value=tip,
                variable=self.method_type,
                command=self.update_method_type,
                style="Custom.TRadiobutton"
            )
            button.grid(row = 1, column = i, padx=10, pady=10)

    def update_str_type(self):
        structure_type = self.str_type.get()
        if structure_type:
            self.structure_type_entry = True
            print('structure_type',structure_type)
            return (structure_type)

    def update_method_type(self):
        method_type = self.method_type.get()
        if method_type:
            self.method_type_entry = True
            print('method_type',method_type)
            return (method_type)

####################################################################################
######## Reading geometry starts here :
####################################################################################

class Read_Geo:
    def __init__(self, file_path):
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
                if element in at_list_bold:
                    atomic_number = at_list_bold.index(element) + 1  # Atomic number = index + 1
                    self.symatno.append(float(atomic_number))
                else:
                    raise ValueError(f"Element {element} is not spelled Correctly.\n"
                                     " or maybe it's above 88 elements of the periodic table \n"
                                     "we only consider fistr 88 elements of the periodic table")


            print("self.symat, self.coordx, self.coordy, self.coordz, self.symatno",self.symat, self.coordx, self.coordy, self.coordz, self.symatno)

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
            global geometry_inserted
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
                if element in at_list_bold:
                    atomic_number = at_list_bold.index(element) + 1  # Atomic number = index + 1
                    self.symatno.append(float(atomic_number))
                else:
                    raise ValueError(f"Element {element} is not spelled Correctly.\n"
                                     " or maybe it's above 88 elements of the periodic table \n"
                                     "we only consider fistr 88 elements of the periodic table")


            print("self.symat, self.coordx, self.coordy, self.coordz, self.symatno",self.symat, self.coordx, self.coordy, self.coordz, self.symatno)
            geometry_inserted = True

            return self.symat, self.coordx, self.coordy, self.coordz, self.symatno


        button = ttk.Button(frame1, text = "Enter", command = create_pan)
        button.grid(row = 0, column = 3, padx = 10, pady = 10, sticky = tk.W )

        button = ttk.Button(frame3, text = "Insert", command = fetch_data)
        button.grid(row = 0, column = 3, padx = 10, pady = 10, sticky = tk.W )
        button = ttk.Button(frame3, text = "Close", command = geo_window.destroy)
        button.grid(row = 1, column = 3, padx = 10, pady = 10, sticky = tk.W )

    def display_geometry(self):
        if geometry_inserted == True:
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


###################################################################################
########## orbital section starts here :
###################################################################################

class Orb_Input:
    def __init__(self,num_orbital, root, keywd_button):
        self.num_orbital = int(num_orbital)

        self.orbital_frame = None
        self.keywd_button = keywd_button
        self.atm_entry = []
        self.typ_entry = []
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
        global type_orb_count
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
        norbsym = ()
        i = 0
        sig_type = 0
        px_type = 0
        py_type = 0
        pz_type = 0
        for data in self.orbital_data:
            i=i+1
            if 's' in data["orbital_type"].lower():
                sig_type = sig_type+1
            elif 'px' in data["orbital_type"].lower():
                px_type = px_type + 1
            elif 'py' in data["orbital_type"].lower():
                py_type = py_type + 1
            elif 'pz' in data["orbital_type"].lower():
                pz_type = pz_type + 1
            elif 'pi1' in data["orbital_type"].lower():
                px_type = px_type + 1
            elif 'pi2' in data["orbital_type"].lower():
                py_type = py_type + 1
            elif 'pi3' in data["orbital_type"].lower():
                pz_type = pz_type + 1
            else:
                messagebox.showerror("Input Error",f"The type of the orbital {i} is unknown. \n"
                                     "Please put s, px, py or pz otherwise put sig or sigma,\n"
                                     " pi1, pi2, pi3 if the direction of pi orbs are different"
                        )
        norbsym = (px_type, py_type, pz_type, sig_type)
        print('norbsym', norbsym)
        print('sig_type, px_type, py_type, pz_type',sig_type, px_type, py_type, pz_type )
        # Check how many types are non-zero
        type_orb_count = sum(1 for count in [sig_type, px_type, py_type, pz_type] if count > 0)

        atoset = self.create_matrix() 
        "atoset is matrix where each row represents an atom according to the geometry, if the atom ia an active atom,"
        "the corresponding column get '1' otherwise '0'. and the next columns contain the corresponding active orbital"
        "numbers associated with that atom."
        print('atoset',atoset)
        self.keywd_button.config(state=tk.NORMAL)  

    def create_matrix(self):
        # Step 1: Gather atom numbers and their orbitals
        print('self.orbital_data',self.orbital_data)
        atom_to_orbitals = {}
        for orbital_number, entry in enumerate(self.orbital_data, start=1):
            atom_number = entry["atom_number"] 
            if atom_number not in atom_to_orbitals:
                atom_to_orbitals[atom_number] = []
            atom_to_orbitals[atom_number].append(orbital_number)
        print('atom_to_orbitals',atom_to_orbitals)

        # Step 2: Determine the maximum atom number
        max_atom_number = max(atom_to_orbitals.keys())
        max_orbitals = max(len(orbitals) for orbitals in atom_to_orbitals.values())

        # Step 3: Initialize the matrix
        matrix = np.zeros((max_atom_number, max_orbitals+1), dtype=int)

        for atom_number, orbitals in atom_to_orbitals.items():
            matrix[atom_number - 1, 0] = 1  # Mark presence
            for col_index, orbital in enumerate(orbitals):
                matrix[atom_number - 1, col_index + 1] = orbital + num_iao # num_iao = number of inactive orbitals

        return matrix


class Keywd_Input:
    def __init__(self, root, multiplicity, type_orb_count, Run_Button):
        self.root = root
        self.multiplicity = multiplicity
        self.type_orb_count = type_orb_count
        self.keywd_window = None
        self.priority_window = None
        self.priority_str_window = None
        self.prio_bond_window = None
        self.prio_rads_window = None
        self.frame = None
        self.set_type = tk.StringVar(value="Single Set")
        self.cheminst_type = tk.StringVar(value="Symmetry")
        self.set_order_type = tk.StringVar(value="Qulity-Arrange")
        self.ovlp_type = tk.StringVar(value="No")
        self.IAB_type = tk.StringVar(value="1")
        self.NAB_type = tk.StringVar(value="2")
        self.SBB_type = tk.StringVar(value="3")
        self.PDB_type = tk.StringVar(value="None")
        self.PDR_type = tk.StringVar(value="None")
        self.set_type_entry = False
        self.cheminst_type_entry = False
        self.set_order_type_entry = False
        self.prio_structure_entries = []  
        self.prio_bond_entries = []  
        self.prio_rads_entries = []  
        self.PDR_buttons = []
        self.SBB_buttons = []
        Pmw.initialise(self.root)
        self.balloon = Pmw.Balloon(self.root)
        mout_number=1

    def create_keywd_pane(self):
        if self.keywd_window is None:
            self.keywd_window = tk.Toplevel(self.root, padx = 10, pady = 10)
            self.keywd_window.title("Spatial Keyword Inputs")
            self.keywd_window.geometry("750x660")
            self.keywd_window.configure( background = "lightblue")

            self.frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.frame.grid (row = 0, column = 0, sticky = tk.W)
            self.frame1 = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.frame1.grid (row = 1, column = 0, sticky = tk.W)
            self.frame2 = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.frame2.grid (row = 5, column = 0, sticky = tk.W)
            self.frame3 = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.frame3.grid (row = 9, column = 0, sticky = tk.W)
            self.frame4 = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.frame4.grid (row = 3, column = 0, sticky = tk.W)
            self.frame5 = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.frame5.grid (row = 4, column = 0, sticky = tk.W)
            self.frame6 = ttk.Frame(self.keywd_window, width = 30,  style="Colour_Frame.TFrame", padding = 10)
            self.frame6.grid(row=6, column=0, sticky="nsew", columnspan = 2)
            self.frame7 = ttk.Frame(self.keywd_window, width = 30,  style="Colour_Frame.TFrame", padding = 10)
            self.frame7.grid(row=7, column=0, sticky="nsew", columnspan = 2)
            self.frame9 = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.frame9.grid (row = 10, column = 0)

            set_label1 = ttk.Label(self.frame, text = "Type of Output Set:", style = "Colour_Label.TLabel")
            set_label1.grid(row = 0, column = 0, sticky = tk.W, padx = 10, pady = 10)
            set_label2 = ttk.Label(self.frame, text = "Symmetric group arrangement:", style = "Colour_Label.TLabel")
            set_label2.grid(row = 1, column = 0, sticky = tk.W, padx = 10, pady = 10)
            set_label3 = ttk.Label(self.frame1, text = "Type of Cheminst Set:", style = "Colour_Label.TLabel")
            set_label3.grid(row = 0, column = 0, padx = 10, pady = 10)

            self.PDB_button = ttk.Button(self.frame7, text = " Insert PDB ", command = self.Insert_PDB)
            self.PDB_button.grid(row = 0, column = 0, padx = 10, sticky = tk.W )
            self.PDB_button.config(state=tk.DISABLED)  # Initially disable the button
            self.PDR_button = ttk.Button(self.frame7, text = " Insert PDR ", command = self.Insert_PDR)
            self.PDR_button.grid(row = 0, column = 1, padx = 10, sticky = tk.W )
            self.PDR_button.config(state=tk.DISABLED)  # Initially disable the button
            insert_button = ttk.Button(self.frame9, text = "Insert", command = self.get_keywds)
            insert_button.grid(row = 0, column = 0, pady = 10, columnspan=2, sticky = tk.W )
            close_button = ttk.Button(self.frame9, text = "Close", command = self.keywd_window.destroy)
            close_button.grid(row = 1, column = 0, pady = 10, columnspan=2, sticky = tk.W )

            Cheminst_type = ["Symmetry", "Asymmetry", "Checksymm"]
            for i, chem in enumerate(Cheminst_type, start=1):
                button1 = ttk.Radiobutton(
                    self.frame,
                    text=chem,
                    value=chem,
                    variable=self.cheminst_type,
                    command=self.cheminst_type_read,
                    style="Custom.TRadiobutton",
                )
                button1.grid(row=0, column=i, padx=10, pady=10)

            set_order = ["Qulity-Arrange","Big-to-Small","Small-to-Big"]
            for i, seto in enumerate(set_order, start=1):
                button3 = ttk.Radiobutton(
                    self.frame,
                    text=seto,
                    value=seto,
                    variable=self.set_order_type,
                    command=self.set_order_type_read,
                    style="Custom.TRadiobutton",
                )
                button3.grid(row=1, column=i, padx=10, pady=10)


            sets_type = ["Single Set", "All Best Sets", "All Sets", "Eq Bond"]
            for i, sets in enumerate(sets_type, start=1):
                button2 = ttk.Radiobutton(
                    self.frame1,
                    text=sets,
                    value=sets,
                    variable=self.set_type,
                    command=self.set_type_read,
                    style="Custom.TRadiobutton",
                )
                button2.grid(row=0, column=i, padx=10, pady=10)


            prio_struc_button = ttk.Button(self.frame3, text = "Insert Priorities Structures", command = self.Insert_Priority_Str )
            prio_struc_button.grid(row = 0, column = 0, padx = 10, pady = 5)
            self.balloon.bind(prio_struc_button, "If you want some structures must be present in the set \n"
                                                 "you can insert them by clicking this Button;\n" 
                                                 "The provided structures must be linearly independent \n"
                                                 "otherwise they will not be present in the set together")
            
            mout_label = ttk.Label(self.frame4, text = "Maximum Number of Output Files:", style = "Colour_Label.TLabel")
            mout_label.grid(row=1,column = 0, padx=10, pady=10, columnspan = 2)
            mout_entry = ttk.Entry(self.frame4, width = 10)
            mout_entry.grid(row = 1, column= 2, padx = 10, pady = 10)
            mout_entry.insert(0, "1")
            self.balloon.bind(mout_label,"For a large Active Space, the total number of sets can\n"
                                         " reach millions or more. Each output file can contain\n" 
                                         " up to 75,000 sets. By default, the number of output files\n"
                                         " is set to one. If additional output files are required, please\n"
                                         " specify the desired number in the entry box and press 'Enter'.")
            
            ovlp_label = ttk.Label(self.frame5, text = "Estimate Overlap", style = "Colour_Label.TLabel")
            ovlp_label.grid(row=0,column = 0, padx=10, pady=10, columnspan = 2)
            self.balloon.bind(ovlp_label, "To estimate overlap among the structures \n"
                                          "in each set click yes. The default is set to 'No'")
            ovlp_types = ["Yes", "No"]
            for i, ovlp in enumerate(ovlp_types, start=2):
                ovlp_button = ttk.Radiobutton(
                    self.frame5,
                    text=ovlp,
                    value=ovlp,
                    variable=self.ovlp_type,
                    command=self.get_str_ovlp_keywd,
                    style="Custom.TRadiobutton",
                )
                ovlp_button.grid(row=0, column=i, padx=10, pady=10)

            prio_label = ttk.Button(self.frame2, text = " Decide Priorities of the Chemical Qualities: ", style = "Colour_Label.TLabel")
            prio_label.grid(row = 0, column = 0, padx = 10, pady = 5)
            self.balloon.bind(prio_label, "The sequence of the default priority is IAB > NAB > SBB and \n"
                              "PBU & PRU are not taken in the calculation; To change this default priorities\n"
                              "please click the button and selevct the numbers for each ")

            IAB_label = ttk.Button(self.frame6, text = " IAB  PRIORITY: ", style = "Colour_Label2.TLabel")
            IAB_label.grid(row = 0, column = 0, padx = 10, pady = 5)
            self.balloon.bind(IAB_label, "Intra Atomic Bond Priority. Default is set to '1'")

            IAB_option = ["1", "2", "3", "4", "5", "None"]
            for i, IAB in enumerate(IAB_option, start=1):
                IAB_button = ttk.Radiobutton(
                    self.frame6,
                    text=IAB,
                    value=IAB,
                    variable=self.IAB_type,
                    command=self.get_IAB_Priority,
                    style="Custom.TRadiobutton",
                )
                IAB_button.grid(row=0, column=i, padx=10, pady=10)

            NAB_label = ttk.Button(self.frame6, text = " NAB PRIORITY: ", style = "Colour_Label2.TLabel")
            NAB_label.grid(row = 1, column = 0, padx = 10, pady = 5)
            self.balloon.bind(NAB_label, "Near Atomic Bond Priority. Default is set to '2'")

            NAB_option = ["1", "2", "3", "4", "5", "None"]
            for i, NAB in enumerate(NAB_option, start=1):
                NAB_button = ttk.Radiobutton(
                    self.frame6,
                    text=NAB,
                    value=NAB,
                    variable=self.NAB_type,
                    command=self.get_NAB_Priority,
                    style="Custom.TRadiobutton",
                )
                NAB_button.grid(row=1, column=i, padx=10, pady=10)

            SBB_label = ttk.Button(self.frame6, text = " SBB PRIORITY: ", style = "Colour_Label2.TLabel")
            SBB_label.grid(row = 2, column = 0, padx = 10, pady = 5)
            self.balloon.bind(SBB_label, "Symmetry Breaking Bond Priority. Default is set to '3'")

            SBB_option = ["1", "2", "3", "4", "5", "None"]
            for i, SBB in enumerate(SBB_option, start=1):
                SBB_button = ttk.Radiobutton(
                    self.frame6,
                    text=SBB,
                    value=SBB,
                    variable=self.SBB_type,
                    command=self.get_SBB_Priority,
                    style="Custom.TRadiobutton",
                )
                SBB_button.grid(row=2, column=i, padx=10, pady=10)
                self.SBB_buttons.append(SBB_button)

            PDB_label = ttk.Button(self.frame6, text = " PDB PRIORITY: ", style = "Colour_Label2.TLabel")
            PDB_label.grid(row = 3, column = 0, padx = 10, pady = 5)
            self.balloon.bind(PDB_label, "Pre-defined Bond Priority. Default is set to \n"
                              "'None': Not taken into the calculation")

            PDB_option = ["1", "2", "3", "4", "5", "None"]
            for i, PDB in enumerate(PDB_option, start=1):
                PDB_button = ttk.Radiobutton(
                    self.frame6,
                    text=PDB,
                    value=PDB,
                    variable=self.PDB_type,
                    command=self.get_PDB_Priority,
                    style="Custom.TRadiobutton",
                )
                PDB_button.grid(row=3, column=i, padx=10, pady=10)

            PDR_label = ttk.Button(self.frame6, text = " PDR PRIORITY: ", style = "Colour_Label2.TLabel")
            PDR_label.grid(row = 4, column = 0, padx = 10, pady = 5)
            self.balloon.bind(PDR_label, "Pre-defined Bond Priority. Default is set to \n"
                              "'None': Not taken into the calculation")

            PDR_option = ["1", "2", "3", "4", "5", "None"]
            for i, PDR in enumerate(PDR_option, start=1):
                PDR_button = ttk.Radiobutton(
                    self.frame6,
                    text=PDR,
                    value=PDR,
                    variable=self.PDR_type,
                    command=self.get_PDR_Priority,
                    style="Custom.TRadiobutton",
                )
                PDR_button.grid(row=4, column=i, padx=10, pady=10)

                self.PDR_buttons.append(PDR_button)
            self.check_and_disable_buttons()


                # Call method to check and disable if multiplicity is 1

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


    def set_type_read(self):
        set_type = self.set_type.get()
        if set_type:
            self.set_type_entry = True
            print('set_type',set_type)
            if set_type == 'All Sets' or set_type =='Eq Bond' or set_type =='All Best Sets':

                ovlp_close_button = ttk.Button(self.frame5, text = "Close", command = self.frame5.destroy)
                ovlp_close_button.grid(row = 0, column= 4, padx = 10, pady = 10)
            return (set_type)

    def get_maximum_num_output(self, mout_entry):
        try:
            mout_number= int(mout_entry.get())
            return(mout_number)
        except ValueError:
            # Handle invalid input gracefully
            tk.messagebox.showerror("Invalid Input", "Please enter a valid number.")
            return
        print('mout_number', mout_number)
#        self.frame4.destroy()

    def get_str_ovlp_keywd(self):
        ovlp_cal = self.ovlp_type.get()
        print('ovlp_cal', ovlp_cal)
        return (ovlp_cal)


    def cheminst_type_read(self):
        cheminst_type = self.cheminst_type.get()
        if cheminst_type:
            self.cheminst_type_entry = True
            print('cheminst_type',cheminst_type)
        return (cheminst_type)
        #    return (cheminst_type)
#            if cheminst_type == 'Symmetry':


    def set_order_type_read(self):
        set_order_type = self.set_order_type.get()
        if set_order_type:
            self.set_order_type_entry = True
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
            bond_number = int(prio_bond_entry.get())  # Get the number of structures
        except ValueError:
            # Handle invalid input gracefully
            tk.messagebox.showerror("Invalid Input", "Please enter a valid number.")
            return

        # Clear previous fields
        for widget in frame1.winfo_children():
            widget.destroy()

        for i in range(bond_number):
            prio_bond_label = ttk.Label(frame1, text=f"Pre-defined Bond {i+1}:", style="Colour_Label1.TLabel")
            prio_bond_label.grid(row=i, column=0, padx=10, pady=10)

            prio_bond_entry = ttk.Entry(frame1, width=20)
            prio_bond_entry.grid(row=i, column=1, padx=10, pady=10)

            self.prio_bond_entries.append(prio_bond_entry)  # Save reference

    def get_prio_bond_data(self):
        """
        Retrieve data from the dynamically created structure entry fields.
        """
        data = [entry.get() for entry in self.prio_bond_entries]
        print("Entered prio_bond:", data)
        self.prio_bond_window.destroy()
        self.prio_bond_window = None
        return data

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
        cheminst = self.cheminst_type_read()
        if cheminst == 'Symmetry':
            symm = 1
        if cheminst == 'Asymmetry':
            symm = 0
        if cheminst == 'Checksymm':
            symtype = 'check'


        sotype = self.set_order_type_read()
        if sotype == 'Quality-Arrange':
            set_order = 0
        elif sotype == 'Small-to-Big':
            set_order = 1
        elif sotype == 'Big-to-Small':
            set_order = 2


        settype = self.set_type_read()
        if settype == 'Single Set':
            nset = 0
        elif settype == 'All Best Sets':
            nset = 1
        elif settype == 'All Sets':
            nset = 2
        elif settype == 'Eq Bond':
            nset = 4
        mout = self.get_maximum_num_output()
        ovlp = self.get_str_ovlp_keywd()
        IAB = self.get_IAB_Priority()
        NAB = self.get_NAB_Priority()
        SBB = self.get_SBB_Priority()
        PDB = self.get_PDB_Priority()
        PDR = self.get_PDR_Priority()
        print("cheminst",cheminst)


def finish():
    sys.exit()
        
def create_orb(inputc, root, keywd_button):
    global num_orbital, num_electron, multiplicity
    num_orbital, num_electron, multiplicity = inputc.generate_ctrl_input()  # calling generate_input to get active orbital number
    num_orbital= int(num_orbital)
    num_electron= int(num_electron)
    multiplicity= int(multiplicity)
    inputo=Orb_Input(num_orbital, root, keywd_button)
    return (num_orbital, num_electron, multiplicity)

def create_keywd(root, Run_button):
    global multiplicity, type_orb_count
    input_keywd = Keywd_Input(root, multiplicity, type_orb_count, Run_button)
    input_keywd.create_keywd_pane()

#def Run_fort_subs(root, ):
    

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Input File Creator")
    root.geometry("550x750")
    root.configure(background="lightblue")
    inputc = Ctrl_Input(root)
    inputc.create_ctrl_pans()

    style_colour_frame = ttk.Style()
    style_colour_frame.configure("Colour_Frame.TFrame", background="lightblue")
    frame1 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
    frame1.grid(row=5, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
    frame2 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
    frame2.grid(row=6, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

    Run_button = ttk.Button(frame1, text="RUN", command=lambda:Run_fort_subs(root))
    Run_button.grid(row=1, column=0, padx=5, pady=10)
    Run_button.config(state=tk.DISABLED)  # Initially disable the button
    inputc.Run_button = Run_button

    keywd_button = ttk.Button(frame1, text="KEYWDS", command=lambda:create_keywd(root, Run_button))
    keywd_button.grid(row=0, column=1, padx=5, pady=10)
    keywd_button.config(state=tk.DISABLED)  # Initially disable the button
    inputc.keywd_button = keywd_button

    orbital_button = ttk.Button(frame1, text="ORBITALS", command=lambda:create_orb(inputc, root, keywd_button))
    orbital_button.grid(row=0, column=0, padx=5, pady=10)
    orbital_button.config(state=tk.DISABLED)  # Initially disable the button
    inputc.orbital_button = orbital_button

    close_button = ttk.Button(frame2, text="FINISH", command=finish)
    close_button.grid(row=0, column=1, pady=10)

    root.mainloop()
