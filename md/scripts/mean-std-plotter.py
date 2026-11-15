#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
import seaborn as sns
import tkinter as tk
from tkinter import filedialog, messagebox, colorchooser, simpledialog
from tkinter import ttk
import matplotlib.colors as mcolors
import pandas as pd
import os
import re

# ============================================================================
# PLOTTING STYLE CONFIGURATION
# ============================================================================

# Font settings - Professional and clean
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
mpl.rcParams['font.size'] = 12

# Figure settings
mpl.rcParams['figure.dpi'] = 100
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['savefig.facecolor'] = 'white'
mpl.rcParams['figure.facecolor'] = 'white'

# Axes settings - Clean and professional
mpl.rcParams['axes.linewidth'] = 1.2
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['axes.titleweight'] = 'bold'
mpl.rcParams['axes.facecolor'] = 'white'
mpl.rcParams['axes.edgecolor'] = 'black'
mpl.rcParams['axes.grid'] = True
mpl.rcParams['axes.grid.axis'] = 'y'
mpl.rcParams['axes.axisbelow'] = True

# Grid settings - Subtle and professional
mpl.rcParams['grid.alpha'] = 0.3
mpl.rcParams['grid.linestyle'] = '--'
mpl.rcParams['grid.linewidth'] = 0.6
mpl.rcParams['grid.color'] = 'gray'

# Lines and markers
mpl.rcParams['lines.linewidth'] = 2.0
mpl.rcParams['lines.markersize'] = 6

# Ticks - Professional appearance
mpl.rcParams['xtick.labelsize'] = 11
mpl.rcParams['ytick.labelsize'] = 11
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.major.width'] = 1.2
mpl.rcParams['ytick.major.width'] = 1.2
mpl.rcParams['xtick.minor.width'] = 0.8
mpl.rcParams['ytick.minor.width'] = 0.8
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['ytick.minor.size'] = 3

# Legend - Clean and readable
mpl.rcParams['legend.fontsize'] = 11
mpl.rcParams['legend.framealpha'] = 0.95
mpl.rcParams['legend.edgecolor'] = 'black'
mpl.rcParams['legend.fancybox'] = False
mpl.rcParams['legend.shadow'] = False

# Professional colorblind-friendly palette
PROFESSIONAL_COLORS = [
    '#0173B2',  # Blue
    '#DE8F05',  # Orange
    '#029E73',  # Teal
    '#CC78BC',  # Purple
    '#ECE133',  # Yellow
    '#56B4E9',  # Light blue
    '#D55E00',  # Red-orange
    '#009E73',  # Green
]

class SecondaryStructureDialog:
    def __init__(self, parent):
        self.result = None
        self.dialog = tk.Toplevel(parent)
        self.dialog.title("Secondary Structure Regions")
        self.dialog.geometry("650x600")
        self.dialog.resizable(True, True)
        
        # Instructions
        instructions = ttk.Label(self.dialog, text="Define secondary structure regions for background shading.\n"
                                                   "Format: start-end (e.g., 10-25 for residues 10 to 25)",
                                justify="left", wraplength=550)
        instructions.pack(padx=10, pady=10)
        
        # Frame for alpha helices
        helix_frame = ttk.LabelFrame(self.dialog, text="Alpha Helices")
        helix_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        ttk.Label(helix_frame, text="Enter alpha-helix regions (one per line):").pack(pady=5)
        self.helix_text = tk.Text(helix_frame, height=5, width=50)
        self.helix_text.pack(padx=10, pady=5)
        self.helix_text.insert("1.0", "# Example:\n# 15-30\n# 45-62\n")
        
        # Color picker for helices
        helix_color_frame = ttk.Frame(helix_frame)
        helix_color_frame.pack(pady=5)
        ttk.Label(helix_color_frame, text="Helix Color:").pack(side="left", padx=5)
        self.helix_color = "#FFB6C1"  # Light pink
        self.helix_color_btn = ttk.Button(helix_color_frame, text="Choose Color", 
                                         command=lambda: self.choose_color('helix'))
        self.helix_color_btn.pack(side="left", padx=5)
        self.helix_color_label = ttk.Label(helix_color_frame, text="  ", background=self.helix_color, width=5)
        self.helix_color_label.pack(side="left", padx=5)
        
        # Frame for beta strands
        strand_frame = ttk.LabelFrame(self.dialog, text="Beta Strands")
        strand_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        ttk.Label(strand_frame, text="Enter beta-strand regions (one per line):").pack(pady=5)
        self.strand_text = tk.Text(strand_frame, height=5, width=50)
        self.strand_text.pack(padx=10, pady=5)
        self.strand_text.insert("1.0", "# Example:\n# 5-12\n# 70-78\n")
        
        # Color picker for strands
        strand_color_frame = ttk.Frame(strand_frame)
        strand_color_frame.pack(pady=5)
        ttk.Label(strand_color_frame, text="Strand Color:").pack(side="left", padx=5)
        self.strand_color = "#ADD8E6"  # Light blue
        self.strand_color_btn = ttk.Button(strand_color_frame, text="Choose Color", 
                                          command=lambda: self.choose_color('strand'))
        self.strand_color_btn.pack(side="left", padx=5)
        self.strand_color_label = ttk.Label(strand_color_frame, text="  ", background=self.strand_color, width=5)
        self.strand_color_label.pack(side="left", padx=5)
        
        # Opacity slider
        opacity_frame = ttk.Frame(self.dialog)
        opacity_frame.pack(fill="x", padx=10, pady=5)
        ttk.Label(opacity_frame, text="Background Opacity:").pack(side="left", padx=5)
        self.opacity_var = tk.DoubleVar(value=0.3)
        self.opacity_slider = ttk.Scale(opacity_frame, from_=0.1, to=0.8, variable=self.opacity_var, 
                                       orient="horizontal", length=200)
        self.opacity_slider.pack(side="left", padx=5)
        self.opacity_label = ttk.Label(opacity_frame, text="0.30")
        self.opacity_label.pack(side="left", padx=5)
        self.opacity_slider.config(command=self.update_opacity_label)
        
        # Buttons
        btn_frame = ttk.Frame(self.dialog)
        btn_frame.pack(fill="x", padx=10, pady=20)
        
        ok_btn = ttk.Button(btn_frame, text="OK", command=self.ok, width=15)
        ok_btn.pack(side="left", padx=10)
        
        cancel_btn = ttk.Button(btn_frame, text="Cancel", command=self.cancel, width=15)
        cancel_btn.pack(side="left", padx=10)
        
        clear_btn = ttk.Button(btn_frame, text="Clear All", command=self.clear_all, width=15)
        clear_btn.pack(side="left", padx=10)
        
        self.dialog.transient(parent)
        self.dialog.grab_set()
        parent.wait_window(self.dialog)
    
    def choose_color(self, structure_type):
        color = colorchooser.askcolor(title=f"Choose {structure_type.capitalize()} Color")
        if color[1]:
            if structure_type == 'helix':
                self.helix_color = color[1]
                self.helix_color_label.config(background=self.helix_color)
            else:
                self.strand_color = color[1]
                self.strand_color_label.config(background=self.strand_color)
    
    def update_opacity_label(self, value):
        self.opacity_label.config(text=f"{float(value):.2f}")
    
    def parse_regions(self, text):
        """Parse text input to extract region ranges"""
        regions = []
        invalid_lines = []
        lines = text.strip().split('\n')
        
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            try:
                # Parse format like "10-25"
                if '-' in line:
                    parts = line.split('-')
                    if len(parts) != 2:
                        invalid_lines.append(f"Line {line_num}: '{line}' (too many dashes)")
                        continue
                    start = int(parts[0].strip())
                    end = int(parts[1].strip())
                    if start >= end:
                        invalid_lines.append(f"Line {line_num}: '{line}' (start must be < end)")
                        continue
                    regions.append((start, end))
                else:
                    invalid_lines.append(f"Line {line_num}: '{line}' (missing dash)")
            except ValueError:
                invalid_lines.append(f"Line {line_num}: '{line}' (invalid numbers)")
        
        # Show warning for invalid lines if any
        if invalid_lines:
            error_msg = "Some lines could not be parsed:\n\n" + "\n".join(invalid_lines[:10])
            if len(invalid_lines) > 10:
                error_msg += f"\n... and {len(invalid_lines) - 10} more"
            messagebox.showwarning("Parsing Warnings", error_msg)
        
        return regions
    
    def ok(self):
        helix_text = self.helix_text.get("1.0", "end")
        strand_text = self.strand_text.get("1.0", "end")
        
        helices = self.parse_regions(helix_text)
        strands = self.parse_regions(strand_text)
        
        if not helices and not strands:
            response = messagebox.askyesno("No Regions Defined", 
                                          "No valid regions defined. Plot will have no shading.\n\n"
                                          "Do you want to continue without secondary structure shading?")
            if not response:
                return  # User wants to go back and define regions
        
        # Show confirmation of what was parsed
        confirmation_msg = "Secondary structure regions defined:\n\n"
        if helices:
            helix_ranges = ", ".join([f"{start}-{end}" for start, end in helices])
            confirmation_msg += f"Alpha Helices ({len(helices)}): {helix_ranges}\n"
        else:
            confirmation_msg += "Alpha Helices: None\n"
        
        if strands:
            strand_ranges = ", ".join([f"{start}-{end}" for start, end in strands])
            confirmation_msg += f"Beta Strands ({len(strands)}): {strand_ranges}\n"
        else:
            confirmation_msg += "Beta Strands: None\n"
        
        confirmation_msg += f"\nHelix color: {self.helix_color}\n"
        confirmation_msg += f"Strand color: {self.strand_color}\n"
        confirmation_msg += f"Opacity: {self.opacity_var.get():.2f}\n\n"
        confirmation_msg += "Is this correct?"
        
        confirmed = messagebox.askyesno("Confirm Secondary Structure", confirmation_msg)
        if not confirmed:
            return  # User wants to go back and edit
        
        self.result = {
            'helices': helices,
            'strands': strands,
            'helix_color': self.helix_color,
            'strand_color': self.strand_color,
            'opacity': self.opacity_var.get()
        }
        self.dialog.destroy()
    
    def cancel(self):
        self.result = None
        self.dialog.destroy()
    
    def clear_all(self):
        self.helix_text.delete("1.0", "end")
        self.strand_text.delete("1.0", "end")

class RMSDPlotterGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("RMSD/RMSF Plotter with Secondary Structure")
        self.root.geometry("1000x850")
        
        # Initialize variables
        self.sets = []
        self.current_plot = None  # To keep track of the current plot
        self.secondary_structure = None  # To store secondary structure info
        
        # Create GUI components
        self.create_widgets()
    
    def create_widgets(self):
        # Frame for adding sets
        set_frame = ttk.LabelFrame(self.root, text="Sets")
        set_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Treeview to display sets
        self.tree = ttk.Treeview(set_frame, columns=("Label", "Color", "Line Style", "Files"), show='headings', height=10)
        self.tree.heading("Label", text="Label")
        self.tree.heading("Color", text="Color")
        self.tree.heading("Line Style", text="Line Style")
        self.tree.heading("Files", text="Number of Files")
        self.tree.column("Label", width=150)
        self.tree.column("Color", width=100)
        self.tree.column("Line Style", width=100)
        self.tree.column("Files", width=150)
        self.tree.pack(fill="both", expand=True, side="left", padx=(0,10), pady=10)
        
        # Scrollbar for treeview
        scrollbar = ttk.Scrollbar(set_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscroll=scrollbar.set)
        scrollbar.pack(side="left", fill="y")
        
        # Buttons for set operations
        btn_frame = ttk.Frame(set_frame)
        btn_frame.pack(fill="y", expand=False, side="left")
        
        add_set_btn = ttk.Button(btn_frame, text="Add Set", command=self.add_set)
        add_set_btn.pack(fill="x", pady=5)
        
        remove_set_btn = ttk.Button(btn_frame, text="Remove Selected Set", command=self.remove_set)
        remove_set_btn.pack(fill="x", pady=5)
        
        # Frame for plot parameters
        param_frame = ttk.LabelFrame(self.root, text="Plot Parameters")
        param_frame.pack(fill="x", expand=False, padx=10, pady=5)
        
        # Axis Labels
        ttk.Label(param_frame, text="X-axis Label:").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.x_label_entry = ttk.Entry(param_frame, width=20)
        self.x_label_entry.grid(row=0, column=1, padx=5, pady=5, sticky="w")
        self.x_label_entry.insert(0, "Time (ns)")
        
        ttk.Label(param_frame, text="Y-axis Label:").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        self.y_label_entry = ttk.Entry(param_frame, width=20)
        self.y_label_entry.grid(row=1, column=1, padx=5, pady=5, sticky="w")
        self.y_label_entry.insert(0, "RMSD (nm)")
        
        # Axis Limits
        ttk.Label(param_frame, text="X-axis Min:").grid(row=0, column=2, padx=5, pady=5, sticky="e")
        self.x_min_entry = ttk.Entry(param_frame, width=10)
        self.x_min_entry.grid(row=0, column=3, padx=5, pady=5, sticky="w")
        self.x_min_entry.insert(0, "0")
        
        ttk.Label(param_frame, text="X-axis Max:").grid(row=1, column=2, padx=5, pady=5, sticky="e")
        self.x_max_entry = ttk.Entry(param_frame, width=10)
        self.x_max_entry.grid(row=1, column=3, padx=5, pady=5, sticky="w")
        self.x_max_entry.insert(0, "100")
        
        ttk.Label(param_frame, text="Y-axis Min:").grid(row=0, column=4, padx=5, pady=5, sticky="e")
        self.y_min_entry = ttk.Entry(param_frame, width=10)
        self.y_min_entry.grid(row=0, column=5, padx=5, pady=5, sticky="w")
        self.y_min_entry.insert(0, "0")
        
        ttk.Label(param_frame, text="Y-axis Max:").grid(row=1, column=4, padx=5, pady=5, sticky="e")
        self.y_max_entry = ttk.Entry(param_frame, width=10)
        self.y_max_entry.grid(row=1, column=5, padx=5, pady=5, sticky="w")
        self.y_max_entry.insert(0, "0.7")
        
        # Tick spacing
        ttk.Label(param_frame, text="X Major Tick:").grid(row=0, column=6, padx=5, pady=5, sticky="e")
        self.x_tick_entry = ttk.Entry(param_frame, width=10)
        self.x_tick_entry.grid(row=0, column=7, padx=5, pady=5, sticky="w")
        self.x_tick_entry.insert(0, "20")
        
        ttk.Label(param_frame, text="Y Major Tick:").grid(row=1, column=6, padx=5, pady=5, sticky="e")
        self.y_tick_entry = ttk.Entry(param_frame, width=10)
        self.y_tick_entry.grid(row=1, column=7, padx=5, pady=5, sticky="w")
        self.y_tick_entry.insert(0, "0.1")
        
        # Secondary structure button
        ss_frame = ttk.LabelFrame(self.root, text="Secondary Structure (for RMSF plots)")
        ss_frame.pack(fill="x", expand=False, padx=10, pady=5)
        
        self.ss_enabled = tk.BooleanVar(value=False)
        ss_check = ttk.Checkbutton(ss_frame, text="Enable secondary structure background shading", 
                                   variable=self.ss_enabled)
        ss_check.pack(side="left", padx=10, pady=5)
        
        ss_btn = ttk.Button(ss_frame, text="Define Secondary Structure Regions", 
                           command=self.define_secondary_structure)
        ss_btn.pack(side="left", padx=10, pady=5)
        
        self.ss_status_label = ttk.Label(ss_frame, text="No regions defined", foreground="gray")
        self.ss_status_label.pack(side="left", padx=10, pady=5)
        
        # Frame for actions
        action_frame = ttk.Frame(self.root)
        action_frame.pack(fill="x", expand=False, padx=10, pady=10)
        
        plot_btn = ttk.Button(action_frame, text="Generate Plot", command=self.generate_plot)
        plot_btn.pack(side="left", padx=5)
        
        save_btn = ttk.Button(action_frame, text="Save Plot", command=self.save_plot)
        save_btn.pack(side="left", padx=5)
        
        clear_btn = ttk.Button(action_frame, text="Clear All Sets", command=self.clear_sets)
        clear_btn.pack(side="left", padx=5)
        
        exit_btn = ttk.Button(action_frame, text="Exit", command=self.root.quit)
        exit_btn.pack(side="right", padx=5)
    
    def define_secondary_structure(self):
        dialog = SecondaryStructureDialog(self.root)
        if dialog.result:
            self.secondary_structure = dialog.result
            n_helices = len(dialog.result['helices'])
            n_strands = len(dialog.result['strands'])
            status_text = f"{n_helices} helix region(s), {n_strands} strand region(s)"
            self.ss_status_label.config(text=status_text, foreground="green")
        else:
            self.ss_status_label.config(text="No regions defined", foreground="gray")
    
    def add_set(self):
        # Open file dialog to select multiple .xvg or .csv files
        files = filedialog.askopenfilenames(title="Select .xvg or .csv Files", 
                                          filetypes=[("XVG files", "*.xvg"), ("CSV files", "*.csv"), ("All files", "*.*")])
        if not files:
            return
        
        # Prompt for label
        label = simpledialog.askstring("Input", "Enter a label for this set:", parent=self.root)
        if not label:
            messagebox.showerror("Error", "Label cannot be empty.")
            return
        
        # Choose color - allow hex input or color picker
        color_input = simpledialog.askstring("Color Selection", 
            "Enter a hex color code (e.g., #FF5733) or leave blank to use color picker:",
            parent=self.root)
        
        if color_input is not None and color_input.strip():
            # Validate hex code
            hex_color = color_input.strip()
            if not hex_color.startswith('#'):
                hex_color = '#' + hex_color
            # Check if valid hex color
            if re.match(r'^#[0-9A-Fa-f]{6}$', hex_color):
                color = (None, hex_color)
            else:
                messagebox.showerror("Error", "Invalid hex color code. Please use format #RRGGBB.")
                return
        else:
            # Use color picker
            color = colorchooser.askcolor(title="Choose Color for Set")
            if color[1] is None:
                messagebox.showerror("Error", "No color selected.")
                return
        
        # Choose line style
        line_styles = ['solid', 'dashed', 'dashdot', 'dotted']
        line_style = simpledialog.askstring("Line Style", 
                                           f"Enter line style for '{label}' set.\nOptions: {', '.join(line_styles)}",
                                           parent=self.root)
        if line_style is None:
            messagebox.showerror("Error", "No line style selected.")
            return
        line_style = line_style.lower()
        if line_style not in line_styles:
            messagebox.showerror("Error", f"Invalid line style selected. Choose from: {', '.join(line_styles)}.")
            return
        
        # Add set to the list
        self.sets.append({
            "label": label,
            "color": color[1],
            "line_style": line_style,
            "files": list(files)
        })
        
        # Insert into treeview
        self.tree.insert("", "end", values=(label, color[1], line_style, len(files)))
    
    def remove_set(self):
        selected_item = self.tree.selection()
        if not selected_item:
            messagebox.showerror("Error", "No set selected.")
            return
        
        index = self.tree.index(selected_item[0])
        del self.sets[index]
        self.tree.delete(selected_item)
    
    def clear_sets(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.sets.clear()
        self.current_plot = None  # Reset the current plot
    
    def read_xvg(self, file_path):
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
            
            # Extract data lines (ignore lines starting with '#' or '@')
            data_lines = [line for line in lines if not line.startswith(('#', '@'))]
            data = np.array([list(map(float, line.split())) for line in data_lines])
            
            if data.shape[1] < 2:
                raise ValueError("File does not contain at least two columns.")
            
            return data
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read {file_path}.\n{e}")
            return None
    
    def read_csv(self, file_path):
        try:
            # Prompt user if CSV has headers
            has_header = messagebox.askyesno("CSV Header", f"Does the file '{os.path.basename(file_path)}' have a header row?")
            delimiter = simpledialog.askstring("Delimiter", f"Enter the delimiter used in '{os.path.basename(file_path)}' (default is comma):", parent=self.root)
            delimiter = delimiter if delimiter else ","
            
            df = pd.read_csv(file_path, delimiter=delimiter, header=0 if has_header else None)
            
            if df.shape[1] < 2:
                raise ValueError("CSV file does not contain at least two columns.")
            
            data = df.values
            return data
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read {file_path}.\n{e}")
            return None
    
    def read_file(self, file_path):
        _, ext = os.path.splitext(file_path)
        if ext.lower() == ".xvg":
            return self.read_xvg(file_path)
        elif ext.lower() == ".csv":
            return self.read_csv(file_path)
        else:
            messagebox.showerror("Error", f"Unsupported file format: {ext}")
            return None
    
    def process_set(self, set_data):
        data_list = []
        x = None
        for file in set_data["files"]:
            data = self.read_file(file)
            if data is None:
                return None, None, None
            if x is None:
                x = data[:, 0]
            else:
                if not np.array_equal(x, data[:, 0]):
                    messagebox.showerror("Error", f"X values in '{os.path.basename(file)}' do not match.")
                    return None, None, None
            data_list.append(data[:, 1])
        
        mean_values = np.mean(data_list, axis=0)
        std_dev_values = np.std(data_list, axis=0)
        return x, mean_values, std_dev_values
    
    def add_secondary_structure_background(self, ax, y_min, y_max):
        """Add shaded background regions for secondary structures"""
        if not self.secondary_structure:
            return
        
        opacity = self.secondary_structure['opacity']
        
        # Add beta strand regions (typically shown in blue)
        for start, end in self.secondary_structure['strands']:
            ax.axvspan(start, end, 
                      color=self.secondary_structure['strand_color'], 
                      alpha=opacity, 
                      zorder=0)
        
        # Add alpha helix regions (typically shown in red/pink)
        for start, end in self.secondary_structure['helices']:
            ax.axvspan(start, end, 
                      color=self.secondary_structure['helix_color'], 
                      alpha=opacity, 
                      zorder=0)
    
    def generate_plot(self):
        if not self.sets:
            messagebox.showerror("Error", "No sets to plot.")
            return
        
        # Get plot parameters
        x_label = self.x_label_entry.get().strip() or "Residue Number"
        y_label = self.y_label_entry.get().strip() or "RMSF (nm)"
        
        try:
            x_min = float(self.x_min_entry.get().strip() or 0)
            x_max = float(self.x_max_entry.get().strip() or 100)
            y_min = float(self.y_min_entry.get().strip() or 0)
            y_max = float(self.y_max_entry.get().strip() or 0.7)
            x_tick = float(self.x_tick_entry.get().strip() or 20)
            y_tick = float(self.y_tick_entry.get().strip() or 0.1)
        except ValueError:
            messagebox.showerror("Error", "Invalid axis limits or tick spacing. Please enter numerical values.")
            return
        
        # Create plot with professional styling
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Add secondary structure background if enabled
        if self.ss_enabled.get() and self.secondary_structure:
            self.add_secondary_structure_background(ax, y_min, y_max)
        
        # Plot data with professional styling
        for idx, set_data in enumerate(self.sets):
            x, mean, std = self.process_set(set_data)
            if x is None:
                return  # Error already shown
            
            # Use professional colors if default colors are being used
            color = set_data["color"]
            if color in ['blue', 'red', 'green', 'orange']:  # Default tkinter colors
                color = PROFESSIONAL_COLORS[idx % len(PROFESSIONAL_COLORS)]
            
            # Plot line with confidence interval
            ax.plot(x, mean, label=f"{set_data['label']}", 
                   color=color, linestyle=set_data["line_style"], 
                   linewidth=1.5, alpha=0.85)
            ax.fill_between(x, mean - std, mean + std, 
                           color=color, alpha=0.2)
        
        # Set labels with professional styling
        ax.set_xlabel(x_label, fontweight='bold', fontsize=12)
        ax.set_ylabel(y_label, fontweight='bold', fontsize=12)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        
        # Professional legend positioning
        #ax.legend(loc='lower-center', bbox_to_anchor=(0.5, 1.02), fontsize=11, framealpha=0.95, 
        #         edgecolor='black', fancybox=False)
        n_sets = len(self.sets)
        legend = ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02),
                           ncol=min(n_sets, 4), frameon=False,
                           fontsize=13, handlelength=1.6, columnspacing=1.8,
                           borderaxespad=0)
        for text in legend.get_texts():
            text.set_fontweight('bold')    

        # Clean axes - remove top and right spines for cleaner look
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
        ax.spines['left'].set_linewidth(1.2)
        ax.spines['bottom'].set_linewidth(1.2)
        
        # Set tick spacing with professional styling
        ax.xaxis.set_major_locator(MultipleLocator(x_tick))
        ax.yaxis.set_major_locator(MultipleLocator(y_tick))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        
        # Configure ticks - outward pointing, cleaner appearance
        ax.tick_params(axis='both', which='major', 
                      width=1.4, length=7, 
                      labelsize=12, direction='out')
        ax.tick_params(axis='both', which='minor', 
                      width=1, length=4, 
                      direction='out')
        #for label in ax.get_xticklabels():
         #   label.set_fontweight('bold')
        #for label in ax.get_yticklabels():
         #   label.set_fontweight('bold')
        
        # Add subtle horizontal grid lines
        ax.grid(True, axis='y', alpha=0.3, linestyle='--', linewidth=0.6)
        ax.grid(True, axis='x', alpha=0.3, linestyle='--', linewidth=0.6)
        ax.set_axisbelow(True)
        
        plt.tight_layout()
        self.current_plot = fig  # Store the current figure
    
        plt.show()
        plt.close(fig)  # Close the plot frame after showing
    
    def save_plot(self):
        if self.current_plot is None:
            messagebox.showerror("Error", "No plot to save. Please generate a plot first.")
            return
        
        # Open save dialog
        file_types = [
            ("PNG", "*.png"),
            ("PDF", "*.pdf"),
            ("SVG", "*.svg"),
            ("EPS", "*.eps"),
            ("All Files", "*.*")
        ]
        save_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=file_types, title="Save Plot As")
        if not save_path:
            return  # User cancelled
        
        # Ask for DPI
        try:
            dpi = simpledialog.askinteger("Input", "Enter resolution (DPI) for the saved plot (e.g., 300):", 
                                          parent=self.root, minvalue=72, maxvalue=1200)
            if dpi is None:
                return  # User cancelled
        except ValueError:
            messagebox.showerror("Error", "Invalid DPI value. Using default 300 DPI.")
            dpi = 300
        
        try:
            self.current_plot.savefig(save_path, dpi=dpi, bbox_inches='tight')
            messagebox.showinfo("Success", f"Plot saved successfully at:\n{save_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save plot.\n{e}")
    
    def clear_sets(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.sets.clear()
        self.current_plot = None  # Reset the current plot

def main():
    root = tk.Tk()
    app = RMSDPlotterGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()