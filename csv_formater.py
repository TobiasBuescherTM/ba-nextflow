import tkinter as tk
import pandas as pd 
from tkinter import ttk, filedialog, messagebox
import csv




def open_csv_file (on_first = True):
    global csv_file
    if on_first:
        csv_file = filedialog.askopenfilename(filetypes=[("CSV files","*.csv")])
      
        df = pd.read_csv(csv_file, sep=",",index_col = 0)
        df.replace(0,"-")
        df.replace("0","-")
        df.loc["qc"] = "FALSE"
        df.loc["Group"] = str(1)
        df.loc["TECH_replica"] = "FALSE"
        
        df.to_csv(csv_file, index = True)
        
        
        
        
    if csv_file:
        try:
            tree.delete(*tree.get_children())
            with open(csv_file, "r" , newline = "", encoding = "utf-8") as f:
                reader = csv.reader(f)
                
                transposed_csv = list(zip(*reader))
                headers = transposed_csv[0]
                
                #headers = next(reader)
                
                
                
                
                
                
                tree["columns"] = headers
                tree["show"] = "headings"
                
                
                
                for i, col in enumerate(headers):
                    tree.heading(col, text=col)
                    tree.column(col, width=120, minwidth=80, anchor="w")
                    
                for col in tree["columns"]:
                    tree.column(col, width=120)
                
                
                    
                for row_num, row in enumerate(transposed_csv[1:], start = 1):
                    tree.insert("", "end", values = row)
                    
            status_label.config(text = f"Loaded: {csv_file}")
            
        except Exception as e:
            messagebox.showerror("Error", str(e))
            

def get_transpose_rows():
    selected = tree.selection()
    selected_data = [tree.item(item)['values'][0] for item in selected]
    
    status_label.config(text=f"Selected {len(selected)} rows")
    
    return selected_data
    
    
    
    
def set_qc():
    col_index = get_transpose_rows()
    df = pd.read_csv(csv_file, sep=",",index_col = 0)
    df.loc["qc",col_index] = "TRUE"
    df.to_csv(csv_file, index = True)
    open_csv_file(on_first = False)


replicate_counter = 1    
def set_tech_replica():
    global replicate_counter
    col_index = get_transpose_rows()
    df = pd.read_csv(csv_file, sep=",",index_col = 0)
    df.loc["TECH_replica",col_index] = str(replicate_counter)
    df.to_csv(csv_file, index = True)
    open_csv_file(on_first = False)
    replicate_counter = replicate_counter + 1

group_counter = 2  
def set_group():
    global group_counter
    col_index = get_transpose_rows()
    df = pd.read_csv(csv_file, sep=",",index_col = 0)
    df.loc["Group",col_index] = str(group_counter)
    df.to_csv(csv_file, index = True)
    open_csv_file(on_first = False)
    group_counter = group_counter + 1
    
def reset_qc():
    df = pd.read_csv(csv_file, sep=",",index_col = 0)
    df.loc["qc"] = "FALSE"
    df.to_csv(csv_file, index = True)
    open_csv_file(on_first = False)
            
def reset_tech_replica():
    df = pd.read_csv(csv_file, sep=",",index_col = 0)
    df.loc["TECH_replica"] = "FALSE"
    df.to_csv(csv_file, index = True)
    open_csv_file(on_first = False)  
    
def reset_group():
    df = pd.read_csv(csv_file, sep=",",index_col = 0)
    df.loc["Group"] = str(1)
    df.to_csv(csv_file, index = True)
    open_csv_file(on_first = False)  
    
# Setup GUI
root = tk.Tk()
root.title("CSV Viewer")
root.geometry("800x600")

title_label = tk.Label(root, 
                      text="Warning: All changes will overwrite selected file", 
                      font=("Arial", 18, "bold"),
                      pady=20)
title_label.pack()


#buttons
# Main horizontal container
main_frame = tk.Frame(root)
main_frame.pack(fill='x', padx=20, pady=20)

# Left side: Two vertical rows
left_container = tk.Frame(main_frame)
left_container.pack(side='left', fill='x', expand=True)

# Row 1: Action buttons (horizontal)
action_frame = tk.Frame(left_container)
action_frame.pack(fill='x', pady=(0,5))

tk.Button(action_frame, text="QC", command=set_qc, width=12).pack(side='left', padx=3)
tk.Button(action_frame, text="Replicate", command=set_tech_replica, width=12).pack(side='left', padx=3)
tk.Button(action_frame, text="Group", command=set_group, width=12).pack(side='left', padx=3)

# Row 2: Reset buttons (horizontal, BELOW actions)
reset_frame = tk.Frame(left_container)
reset_frame.pack(fill='x', pady=(0,0))

tk.Button(reset_frame, text="Reset QC", command=reset_qc, width=12).pack(side='left', padx=3)
tk.Button(reset_frame, text="Reset Replicate", command=reset_tech_replica, width=12).pack(side='left', padx=3)
tk.Button(reset_frame, text="Reset Group", command=reset_group, width=12).pack(side='left', padx=3)

# Right side: File buttons (unchanged)
right_frame = tk.Frame(main_frame)
right_frame.pack(side='right', padx=(10,0))

tk.Button(right_frame, text="Open CSV", command=open_csv_file).pack(pady=5)
tk.Button(right_frame, text="Quit", command=root.destroy, bg='red', fg='white').pack(pady=5)


# Create scrollbars + frame
tree_frame = tk.Frame(root)
tree_frame.pack(fill='both', expand=True, padx=10, pady=20)
v_scrollbar = ttk.Scrollbar(tree_frame, orient="vertical")
h_scrollbar = ttk.Scrollbar(tree_frame, orient="horizontal")

tree = ttk.Treeview(tree_frame, 
                   yscrollcommand=v_scrollbar.set, 
                   xscrollcommand=h_scrollbar.set)


def drag_select(event):
    item = tree.identify_row(event.y)
    if item:
        tree.selection_set(item)

def extend_select(event):
    item = tree.identify_row(event.y)
    if item:
        tree.selection_add(item)

tree.bind("<Button-1>", drag_select)
tree.bind("<B1-Motion>", extend_select)

v_scrollbar.config(command=tree.yview)
h_scrollbar.config(command=tree.xview)
v_scrollbar.pack(side="right", fill="y")
h_scrollbar.pack(side="bottom", fill="x")
tree.pack(side="left", fill="both", expand=True)

status_label = tk.Label(root, text="No file loaded")
status_label.pack(pady=5)

root.mainloop()