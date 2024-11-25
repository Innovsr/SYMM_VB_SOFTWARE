################################################################################################################
# python GUI interface to create the input file to run the fortran symmetric structure generation symm_str code#
################################################################################################################

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import re
import os
#import program_main
import subprocess
import sys


class Input_creater:
    def __init__(self, root):
        self.root = root
        self.root.title("Input File Creator")

# initialise disctionary to store control data
        self.ctrl_inputs={}
        self.input_text=''
        self.fname=''
 
        # Main Frame
        self.frame = ttk.Frame(root, padding="10")
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Input Fields
        self.entries = {}
        self.keywords = ["str", "nao", "nae", "nmul", "frgtyp", "chinst"]
        self.create_input_fields()

        # Buttons
        self.create_buttons()

        # Fragment Input Pane
        self.fragment_frame = None
        self.fragment_entries = []
        self.description = ""

        self.orbital_frame = None
        self.orbital_entries = []
        self.description_orb = ""

        self.keywd_frame = None
        self.keywd_entries = []
        self.description_key = ""

    def create_input_fields(self):
        ttk.Label(self.frame, text="Enter Ctrl Keywords").grid(row=1, column=0, columnspan=2, pady=5)
        ttk.Label(self.frame, text="file_name").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.file_name = ttk.Entry(self.frame, width=20)
        self.file_name.grid(row=0, column=1, padx=5, pady=5)
#        self.entries["file_name"]=entry
        for idx, key in enumerate(self.keywords):
            ttk.Label(self.frame, text=key).grid(row=idx + 2, column=0, sticky=tk.W, padx=5, pady=5)
            entry = ttk.Entry(self.frame, width=20)
            entry.grid(row=idx + 2, column=1, padx=5, pady=5)
            self.entries[key] = entry


    def create_buttons(self):
        fragment_button = ttk.Button(self.frame, text="Fragments", command=self.create_fragment_section)
        fragment_button.grid(row=8, column=0, padx=5, pady=10)

        orbital_button = ttk.Button(self.frame, text="Orbitals", command=self.create_orbital_section)
        orbital_button.grid(row=8, column=1, padx=5, pady=10)

        keywd_button = ttk.Button(self.frame, text="spatial keywords", command=self.create_keywd_section)
        keywd_button.grid(row=8, column=2, padx=5, pady=10)

        generate_button = ttk.Button(self.frame, text="Generate Input", command=self.generate_input)
        generate_button.grid(row=12, column=0, padx=5, pady=10, sticky=tk.E)

        reset_button = ttk.Button(self.frame, text="Reset", command=self.reset_fields)
        reset_button.grid(row=12, column=2, padx=5, pady=10, sticky=tk.W)


    def go_to_output(self):
        self.frame.destroy()  # Destroy the output frame
        InputCreator(self.root)  # Switch back to the input frame


    def reset_fields(self):
        for entry in self.entries.values():
            entry.delete(0, tk.END)
        if self.fragment_entries:
            for pane in self.fragment_entries:
                for field in pane:
                    field.delete(0, tk.END)
        if self.orbital_entries:
            for pane in self.orbital_entries:
                for field in pane:
                    field.delete(0, tk.END)

################################################################################
#### fragment section starts here : first ##
################################################################################
    def create_fragment_section(self):
        # Create a new frame for fragment inputs if it doesn't exist
        if self.fragment_frame is None:
            self.fragment_frame = tk.Toplevel(self.root, padx=10, pady=10)
            self.fragment_frame.title("Fragment inputs")
            self.fragment_frame.geometry("500x300")

#        # Clear existing fragment panes
#        for widget in self.fragment_frame.winfo_children():
#            widget.destroy()
        self.fragment_entries = []

        # Input for description of fragments
        ttk.Label(self.fragment_frame, text="Description of Fragments").grid(
            row=0, column=2, columnspan=2, pady=5
        )
        desc_entry = ttk.Entry(self.fragment_frame, width=60)
        desc_entry.grid(row=1, column=2, columnspan=2, pady=5)
        ttk.Button(
            self.fragment_frame, text="Analyze Description", command=lambda: self.analyze_fragments(desc_entry.get())
        ).grid(row=2, column=2, columnspan=2, pady=10)

    def analyze_fragments(self, description):
        self.description = description
        try:
            # Calculate the total number of fragments
            num_fragments = self.calculate_fragments(description)
            if num_fragments <= 0:
                raise ValueError("Number of fragments must be positive.")
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter a valid description of fragments.")
            return

        # Clear existing fragment panes
        for widget in self.fragment_frame.winfo_children():
            widget.destroy()

        # Create input panes for each fragment
        ttk.Label(self.fragment_frame, text=f"Number of Fragments: {num_fragments}").grid(
            row=0, column=0, columnspan=2, pady=5
        )
        for i in range(num_fragments):
            row = i + 1
            ttk.Label(self.fragment_frame, text=f"Fragment {i + 1} Type").grid(row=row, column=0, padx=5, pady=5, sticky=tk.W)
            type_entry = ttk.Entry(self.fragment_frame, width=10)
            type_entry.grid(row=row, column=1, padx=5, pady=5)

            ttk.Label(self.fragment_frame, text="Atom Number").grid(row=row, column=2, padx=5, pady=5, sticky=tk.W)
            number_entry = ttk.Entry(self.fragment_frame, width=10)
            number_entry.grid(row=row, column=3, padx=5, pady=5)

            self.fragment_entries.append((type_entry, number_entry))

    def calculate_fragments(self, description):
        """
        Parse the description and calculate the total number of fragments.
        """
        total_fragments = 0
        groups = description.split()
        for group in groups:
            if '*' in group:
                # Handle compact notation like "3*4"
                match = re.match(r"(\d+)\*(\d+)", group)
                if match:
                    total_fragments += int(match.group(2))
            else:
                # Handle individual groups
                total_fragments += 1
        return total_fragments

###################################################################################
########## orbital section starts here :
###################################################################################

    def create_orbital_section(self):
        # Create a new frame for orbital inputs if it doesn't exist
        if self.orbital_frame is None:
            self.orbital_frame = tk.Toplevel(self.root, padx=10, pady=10)
            self.orbital_frame.title("orbital inputs")
            self.orbital_frame.geometry("600x600")

        # Clear existing fragment panes
        for widget in self.orbital_frame.winfo_children():
            widget.destroy()
        self.orbital_entries = []

        # Input for description of fragments
        ttk.Label(self.orbital_frame, text="Description of orbitals").grid(
            row=0, column=0, columnspan=2, pady=5
        )
        desc_entry = ttk.Entry(self.orbital_frame, width=40)
        desc_entry.grid(row=1, column=0, columnspan=2, pady=5)
        ttk.Button(
            self.orbital_frame, text="Analyze Description", command=lambda: self.analyze_orbital(desc_entry.get())
        ).grid(row=2, column=0, columnspan=2, pady=10)


    def analyze_orbital(self, description_orb):
        self.description_orb = description_orb
        try:
            # Calculate the total number of orbitals
            num_orbital = self.calculate_orbital(description_orb)
            if num_orbital <= 0:
                raise ValueError("Number of orbitals must be positive.")
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter a valid description of orbitals.")
            return

        # Clear existing fragment panes
        for widget in self.orbital_frame.winfo_children():
            widget.destroy()

        # Create input panes for each fragment
        ttk.Label(self.orbital_frame, text=f"Number of orbitals: {num_orbital}").grid(
            row=0, column=0, columnspan=2, pady=5
        )
        if num_orbital % 2 == 0:
            j1=int(num_orbital/2)
            j2=int(num_orbital/2)
        else:
            j1=int(num_orbital/2)
            j2=int(num_orbital/2)+1

        j=0
        for i in range(j2):
            j = j+1
            row = i+1
            ttk.Label(self.orbital_frame, text=f"orbital {j } ").grid(row=row, column=1, padx=5, pady=5, sticky=tk.W)
            number_entry = ttk.Entry(self.orbital_frame, width=10)
            number_entry.grid(row=row, column=2, padx=5, pady=5)

            self.orbital_entries.append(number_entry)
        for i in range(j1):
            j = j+1
            row = i+1
            ttk.Label(self.orbital_frame, text=f"orbital {j } ").grid(row=row, column=3, padx=5, pady=5, sticky=tk.W)
            number_entry = ttk.Entry(self.orbital_frame, width=10)
            number_entry.grid(row=row, column=4, padx=5, pady=5)

            self.orbital_entries.append(number_entry)


    def calculate_orbital(self, description_orb):
        """
        Parse the description and calculate the total number of fragments.
        """
        total_orbital = 0
        groups = description_orb.split()
        for group in groups:
            if '*' in group:
                # Handle compact notation like "3*4"
                match = re.match(r"(\d+)\*(\d+)", group)
                if match:
                    total_orbital += int(match.group(2))
            else:
                # Handle individual groups
                total_orbital += 1
        return total_orbital
################################################################################
#### Spatial keyword section starts here : first ##
################################################################################
    def create_keywd_section(self):
        # Create a new frame for fragment inputs if it doesn't exist
        if self.keywd_frame is None:
            self.keywd_frame = tk.Toplevel(self.root, padx=10, pady=10)
            self.keywd_frame.title("Spatial Keywords")
            self.keywd_frame.geometry("400x300")

        # Input for description of fragments
        ttk.Label(self.keywd_frame, text="Number of Spatial Inputs").grid(row=0, column=0, columnspan=2, pady=5)
        desc_entry = ttk.Entry(self.keywd_frame, width=40)
        desc_entry.grid(row=1, column=0, columnspan=2, pady=5)
        ttk.Button(
            self.keywd_frame, text="Create Inputs", command=lambda: self.analyze_keywds(desc_entry.get())
        ).grid(row=2, column=0, columnspan=2, pady=10)

    def analyze_keywds(self, description):
        self.description_key = int(description)
        try:
            # Calculate the total number of fragments
            num_keywds = self.description_key
            if num_keywds <= 0:
                raise ValueError("Number of keywords must be positive.")
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter a valid number of keywords.")
            return

        # Clear existing key word panes
        for widget in self.keywd_frame.winfo_children():
            widget.destroy()

        # Create input panes for each fragment
        ttk.Label(self.keywd_frame, text=f"Number of keywords: {num_keywds}").grid(
            row=0, column=0, columnspan=2, pady=5
        )
        for i in range(num_keywds):
            row = i + 1
            ttk.Label(self.keywd_frame, text=f"keywords {i + 1} Type").grid(row=row, column=0, padx=5, pady=5, sticky=tk.W)
            type_entry = ttk.Entry(self.keywd_frame, width=10)
            type_entry.grid(row=row, column=1, padx=5, pady=5)

            ttk.Label(self.keywd_frame, text="keywords").grid(row=row, column=2, padx=5, pady=5, sticky=tk.W)
            number_entry = ttk.Entry(self.keywd_frame, width=10)
            number_entry.grid(row=row, column=3, padx=5, pady=5)

            self.keywd_entries.append((type_entry, number_entry))

    def show_result(self, input_text):
        result_window = tk.Toplevel(self.root)
        result_window.title("Generated Input")
        result_window.geometry("400x300")

        text_box = tk.Text(result_window, wrap=tk.WORD)
        text_box.insert(tk.END, input_text)
        text_box.config(state=tk.DISABLED)
        text_box.pack(expand=True, fill=tk.BOTH, padx=10, pady=10)

        ttk.Button(result_window, text="Close", command=result_window.destroy).pack(pady=5)

    def generate_input(self):
        # Collect standard inputs
        self.fname=self.file_name.get()
        for key, entry in self.entries.items():
            try:
                value = entry.get()
                if value:
                    self.ctrl_inputs[key]=value
            except AttributeError:
                raise TypeError(f"Expected tk.Entry, but got {type(entry)} for key {key}")
            print('inputs',self.ctrl_inputs)

        # Collect fragment inputs if available
        fragment_inputs = [self.description]
        for pane in self.fragment_entries:
            fragment_type = pane[0].get()
            fragment_number = pane[1].get()
            if fragment_type and fragment_number:
                fragment_inputs.append((fragment_type, fragment_number))
            else:
                raise ValueError(f"Invalid pane in fragment_entries: {pane}")

        # collect orbital inputs if available
        orbital_inputs = [self.description_orb]
        for pane in self.orbital_entries:
            orbital_number = pane.get()
            if orbital_number:
                orbital_inputs.append(orbital_number)
            else:
                raise ValueError(f"Invalid pane in orbital_entries: {pane}")

        # collect spatial keyword inputs if available
        keywd_inputs = [self.description_key]
        for pane in self.keywd_entries:
            keywd_type = pane[0].get()
            keywd_number = pane[1].get()
            if keywd_type and keywd_number:
                keywd_inputs.append((keywd_type, keywd_number))
            else:
                raise ValueError(f"Invalid pane in keyword_entries: {pane}")

        # Generate input string
        self.input_text = f"{self.fname}\n"
        self.input_text += "$ctrl\n"
        self.input_text += " ".join(f"{key}={value}" for key, value in self.ctrl_inputs.items() if key != "fragments" and key != "orbital")
        self.input_text += "\n$end"
        self.input_text += "\n$frag\n"
        self.input_text += f"{self.description}\n"
        for frag in fragment_inputs[1:]:
            self.input_text += f"{frag[0]} {frag[1]}\n"
        self.input_text += "$end"
        self.input_text += "\n$orb\n"
        self.input_text += f"{self.description_orb}\n"
        for orb in orbital_inputs[1:]:
            self.input_text += f"{orb}\n"
        self.input_text += "$end\n"
        for keywd in keywd_inputs[1:]:
            self.input_text += f"{keywd[0]} {keywd[1]}\n"

        self.show_result(self.input_text)
    

    def get_data(self):
        input_text1=self.input_text
        return input_text1
    def get_file_name(self):
        print('file_name',self.fname)
        return self.fname

class analyse_inputs:
    def __init__(self, input_data):
        self.input_data=input_data

    def creat_input(self, file_name1):
        # Write the input text to a file
        self.file_name=''.join([file_name1, ".xmi"])
        with open(self.file_name, 'w') as f:
            f.write(self.input_data)
    def run_Fortran(self):
        #Calling the Fortran executable using subprocess
        process = subprocess.Popen(['./my_program',self.file_name], 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
        
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            print(f"Error running Fortran program: {stderr.decode()}")
        else:
            print(f"Fortran program executed successfully.\nOutput:\n{stdout.decode()}")

class Output:
    def __init__(self, root):
        # Initialize the root window
        self.root = root
        self.root.title("Output Value")

        # Create a frame for the widgets
        self.frame = ttk.Frame(root, padding="10")
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Create a frame for the text widget and scrollbar
        text_frame = ttk.Frame(self.frame)
        text_frame.grid(row=1, column=0, pady=(0, 10), sticky=(tk.W, tk.E, tk.N, tk.S))

        # Add a vertical scrollbar
        self.scrollbar = ttk.Scrollbar(text_frame, orient=tk.VERTICAL)
        self.scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))

        # Add a Text widget to display the file content
        self.output_text = tk.Text(
            text_frame, wrap=tk.WORD, height=25, width=80, yscrollcommand=self.scrollbar.set
        )
        self.output_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Configure the scrollbar to work with the Text widget
        self.scrollbar.config(command=self.output_text.yview)

    def display_file_content(self, output_file):
        """Displays the content of the specified file in the Text widget."""
        with open(output_file, 'r') as f:
            content = f.read()

        # Display the content in the Text widget
        self.output_text.delete(1.0, tk.END)  # Clear previous content
        self.output_text.insert(tk.END, content)

def open_second_root():
    # getting inputs and trasfer for analyse and creating input file
    input_data=inputc.get_data()
    file_name=inputc.get_file_name()
    result=analyse_inputs(input_data)
    result.creat_input(file_name)
    result.run_Fortran()     # instance to run fortran code

    # opening new Tkinter window and showing the results
    root1 = tk.Tk()
    output = Output(root1)
    output_file = "structure_set_1.dat"
    output.display_file_content(output_file)
    root1.mainloop()


def finish():
    sys.exit()


if __name__ == "__main__":
    root = tk.Tk()
    inputc = Input_creater(root)
    run_button = ttk.Button(root, text="RUN", command=open_second_root)
    run_button.grid(row=1, column=0, pady=10)
    close_button = ttk.Button(root, text="FINISH", command=finish)
    close_button.grid(row=2, column=0, pady=10)

    root.mainloop()
#if __name__ == "__main__":
