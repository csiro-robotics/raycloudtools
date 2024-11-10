import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import subprocess
import os

class RayCloudToolsGUI:
    def __init__(self, master):
        self.master = master
        master.title("RayCloudTools GUI for Apptainer")
        master.geometry("800x900")

        self.notebook = ttk.Notebook(master)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.create_rayimport_tab()
        self.create_rayextract_tab()

        tk.Button(master, text="Run Command", command=self.run_command).pack(pady=10)

        self.output_text = tk.Text(master, height=10, width=80)
        self.output_text.pack(pady=10)

    def create_rayimport_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="rayimport")

        tk.Label(tab, text="Point cloud file:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
        self.pointcloud_file = tk.Entry(tab, width=50)
        self.pointcloud_file.grid(row=0, column=1, padx=5, pady=5)
        tk.Button(tab, text="Browse", command=lambda: self.browse_file(self.pointcloud_file)).grid(row=0, column=2, padx=5, pady=5)

        tk.Label(tab, text="Trajectory file:").grid(row=1, column=0, sticky="e", padx=5, pady=5)
        self.trajectory_file = tk.Entry(tab, width=50)
        self.trajectory_file.grid(row=1, column=1, padx=5, pady=5)
        tk.Button(tab, text="Browse", command=lambda: self.browse_file(self.trajectory_file)).grid(row=1, column=2, padx=5, pady=5)

        tk.Label(tab, text="Transform file:").grid(row=2, column=0, sticky="e", padx=5, pady=5)
        self.transform_file = tk.Entry(tab, width=50)
        self.transform_file.grid(row=2, column=1, padx=5, pady=5)
        tk.Button(tab, text="Browse", command=lambda: self.browse_file(self.transform_file)).grid(row=2, column=2, padx=5, pady=5)

        tk.Label(tab, text="Constant position:").grid(row=3, column=0, sticky="e", padx=5, pady=5)
        self.constant_position = tk.Entry(tab, width=50)
        self.constant_position.grid(row=3, column=1, padx=5, pady=5)

        tk.Label(tab, text="Constant ray:").grid(row=4, column=0, sticky="e", padx=5, pady=5)
        self.constant_ray = tk.Entry(tab, width=50)
        self.constant_ray.grid(row=4, column=1, padx=5, pady=5)

        tk.Label(tab, text="Max intensity:").grid(row=5, column=0, sticky="e", padx=5, pady=5)
        self.max_intensity = tk.Entry(tab, width=50)
        self.max_intensity.grid(row=5, column=1, padx=5, pady=5)

        self.remove_start_pos = tk.BooleanVar()
        tk.Checkbutton(tab, text="Remove Start Position", variable=self.remove_start_pos).grid(row=6, column=1, sticky="w", padx=5, pady=5)

    def create_rayextract_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="rayextract")

        self.rayextract_notebook = ttk.Notebook(tab)
        self.rayextract_notebook.pack(fill=tk.BOTH, expand=True)

        self.create_terrain_tab()
        self.create_trunks_tab()
        self.create_forest_tab()
        self.create_trees_tab()
        self.create_leaves_tab()
        self.create_grid_tab()

    def create_terrain_tab(self):
        tab = ttk.Frame(self.rayextract_notebook)
        self.rayextract_notebook.add(tab, text="Terrain")

        tk.Label(tab, text="Cloud file:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
        self.terrain_cloud_file = tk.Entry(tab, width=50)
        self.terrain_cloud_file.grid(row=0, column=1, padx=5, pady=5)
        tk.Button(tab, text="Browse", command=lambda: self.browse_file(self.terrain_cloud_file)).grid(row=0, column=2, padx=5, pady=5)

        tk.Label(tab, text="Gradient:").grid(row=1, column=0, sticky="e", padx=5, pady=5)
        self.terrain_gradient = tk.Entry(tab, width=10)
        self.terrain_gradient.grid(row=1, column=1, sticky="w", padx=5, pady=5)

    def create_trunks_tab(self):
        tab = ttk.Frame(self.rayextract_notebook)
        self.rayextract_notebook.add(tab, text="Trunks")

        tk.Label(tab, text="Cloud file:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
        self.trunks_cloud_file = tk.Entry(tab, width=50)
        self.trunks_cloud_file.grid(row=0, column=1, padx=5, pady=5)
        tk.Button(tab, text="Browse", command=lambda: self.browse_file(self.trunks_cloud_file)).grid(row=0, column=2, padx=5, pady=5)

        self.trunks_exclude_rays = tk.BooleanVar()
        tk.Checkbutton(tab, text="Exclude rays", variable=self.trunks_exclude_rays).grid(row=1, column=1, sticky="w", padx=5, pady=5)

    def create_forest_tab(self):
        tab = ttk.Frame(self.rayextract_notebook)
        self.rayextract_notebook.add(tab, text="Forest")

        tk.Label(tab, text="Cloud file:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
        self.forest_cloud_file = tk.Entry(tab, width=50)
        self.forest_cloud_file.grid(row=0, column=1, padx=5, pady=5)
        tk.Button(tab, text="Browse", command=lambda: self.browse_file(self.forest_cloud_file)).grid(row=0, column=2, padx=5, pady=5)

        tk.Label(tab, text="Ground mesh:").grid(row=1, column=0, sticky="e", padx=5, pady=5)
        self.forest_ground_mesh = tk.Entry(tab, width=50)
        self.forest_ground_mesh.grid(row=1, column=1, padx=5, pady=5)
        tk.Button(tab, text="Browse", command=lambda: self.browse_file(self.forest_ground_mesh)).grid(row=1, column=2, padx=5, pady=5)

        # Add more forest options...

    def create_trees_tab(self):
        tab = ttk.Frame(self.rayextract_notebook)
        self.rayextract_notebook.add(tab, text="Trees")

        tk.Label(tab, text="Cloud file:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
        self.trees_cloud_file = tk.Entry(tab, width=50)
        self.trees_cloud_file.grid(row=0, column=1, padx=5, pady=5)
        tk.Button(tab, text="Browse", command=lambda: self.browse_file(self.trees_cloud_file)).grid(row=0, column=2, padx=5, pady=5)

        tk.Label(tab, text="Ground mesh:").grid(row=1, column=0, sticky="e", padx=5, pady=5)
        self.trees_ground_mesh = tk.Entry(tab, width=50)
        self.trees_ground_mesh.grid(row=1, column=1, padx=5, pady=5)
        tk.Button(tab, text="Browse", command=lambda: self.browse_file(self.trees_ground_mesh)).grid(row=1, column=2, padx=5, pady=5)

        # Add more trees options...

    def create_leaves_tab(self):
        tab = ttk.Frame(self.rayextract_notebook)
        self.rayextract_notebook.add(tab, text="Leaves")

        tk.Label(tab, text="Cloud file:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
        self.leaves_cloud_file = tk.Entry(tab, width=50)
        self.leaves_cloud_file.grid(row=0, column=1, padx=5, pady=5)
        tk.Button(tab, text="Browse", command=lambda: self.browse_file(self.leaves_cloud_file)).grid(row=0, column=2, padx=5, pady=5)

        tk.Label(tab, text="Trees file:").grid(row=1, column=0, sticky="e", padx=5, pady=5)
        self.leaves_trees_file = tk.Entry(tab, width=50)
        self.leaves_trees_file.grid(row=1, column=1, padx=5, pady=5)
        tk.Button(tab, text="Browse", command=lambda: self.browse_file(self.leaves_trees_file)).grid(row=1, column=2, padx=5, pady=5)

        # Add more leaves options...

    def create_grid_tab(self):
        tab = ttk.Frame(self.rayextract_notebook)
        self.rayextract_notebook.add(tab, text="Grid")

        tk.Label(tab, text="Cloud file:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
        self.grid_cloud_file = tk.Entry(tab, width=50)
        self.grid_cloud_file.grid(row=0, column=1, padx=5, pady=5)
        tk.Button(tab, text="Browse", command=lambda: self.browse_file(self.grid_cloud_file)).grid(row=0, column=2, padx=5, pady=5)

        tk.Label(tab, text="Voxel size:").grid(row=1, column=0, sticky="e", padx=5, pady=5)
        self.grid_voxel_size = tk.Entry(tab, width=10)
        self.grid_voxel_size.grid(row=1, column=1, sticky="w", padx=5, pady=5)

        # Add more grid options...

    def browse_file(self, entry):
        filename = filedialog.askopenfilename()
        if filename:
            filename = os.path.abspath(filename)
            entry.delete(0, tk.END)
            entry.insert(0, filename)

    def run_command(self):
        command = ["apptainer", "run", "./raycloudtools_private.sif"]
        
        active_tab = self.notebook.tab(self.notebook.select(), "text").lower()
        
        if active_tab == "rayimport":
            command.append("rayimport")
            command.append(self.pointcloud_file.get())
            if self.trajectory_file.get():
                command.append(self.trajectory_file.get())
            elif self.transform_file.get():
                command.extend(["transform", self.transform_file.get()])
            elif self.constant_position.get():
                command.append(self.constant_position.get())
            elif self.constant_ray.get():
                command.extend(["ray", self.constant_ray.get()])
            
            if self.max_intensity.get():
                command.extend(["--max_intensity", self.max_intensity.get()])
            if self.remove_start_pos.get():
                command.append("--remove_start_pos")
        
        elif active_tab == "rayextract":
            command.append("rayextract")
            active_subtab = self.rayextract_notebook.tab(self.rayextract_notebook.select(), "text").lower()
            
            if active_subtab == "terrain":
                command.extend(["terrain", self.terrain_cloud_file.get()])
                if self.terrain_gradient.get():
                    command.extend(["--gradient", self.terrain_gradient.get()])
            elif active_subtab == "trunks":
                command.extend(["trunks", self.trunks_cloud_file.get()])
                if self.trunks_exclude_rays.get():
                    command.append("--exclude_rays")
            elif active_subtab == "forest":
                command.extend(["forest", self.forest_cloud_file.get()])
                if self.forest_ground_mesh.get():
                    command.extend(["--ground", self.forest_ground_mesh.get()])
                # Add more forest options...
            elif active_subtab == "trees":
                command.extend(["trees", self.trees_cloud_file.get(), self.trees_ground_mesh.get()])
                # Add more trees options...
            elif active_subtab == "leaves":
                command.extend(["leaves", self.leaves_cloud_file.get(), self.leaves_trees_file.get()])
                # Add more leaves options...
            elif active_subtab == "grid":
                command.extend(["grid", self.grid_cloud_file.get()])
                if self.grid_voxel_size.get():
                    command.extend(["--voxel_size", self.grid_voxel_size.get()])
                # Add more grid options...

        # Clear previous output
        self.output_text.delete('1.0', tk.END)

        # Display the command
        self.output_text.insert(tk.END, f"Running command:\n{' '.join(command)}\n\n")
        self.output_text.update()

        try:
            # Run the command and capture output in real-time
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            
            for line in process.stdout:
                self.output_text.insert(tk.END, line)
                self.output_text.see(tk.END)
                self.output_text.update()
            
            process.wait()
            
            if process.returncode == 0:
                messagebox.showinfo("Success", "Command completed successfully")
            else:
                messagebox.showerror("Error", f"Command failed with return code {process.returncode}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to run command: {str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    gui = RayCloudToolsGUI(root)
    root.mainloop()
