import py3Dmol
import numpy as np

def addBox(view, vecs, label=True, viewer=None):
    pts = [
           ([0,0,0], vecs[0,:]),
           ([0,0,0], vecs[1,:]),
           ([0,0,0], vecs[2,:]),
           (vecs[0,:], vecs[0,:]+vecs[1,:]),
           (vecs[0,:], vecs[0,:]+vecs[2,:]),
           (vecs[1,:], vecs[0,:]+vecs[1,:]),
           (vecs[1,:], vecs[2,:]+vecs[1,:]),
           (vecs[2,:], vecs[0,:]+vecs[2,:]),
           (vecs[2,:], vecs[1,:]+vecs[2,:]),
           (vecs[0,:]+vecs[1,:], vecs[0,:]+vecs[1,:]+vecs[2,:]),
           (vecs[0,:]+vecs[2,:], vecs[0,:]+vecs[1,:]+vecs[2,:]),
           (vecs[1,:]+vecs[2,:], vecs[0,:]+vecs[1,:]+vecs[2,:]),
          ]

    for i, pt in enumerate(pts):
        pt1, pt2 = pt

        if viewer is None:
            view.addLine({"start": {"x":pt1[0], "y": pt1[1], "z": pt1[2]}, 
                           "end":  {"x":pt2[0], "y": pt2[1], "z": pt2[2]}})
        else:
            view.addLine({"start": {"x":pt1[0], "y": pt1[1], "z": pt1[2]}, 
                           "end":  {"x":pt2[0], "y": pt2[1], "z": pt2[2]}},
                           viewer = viewer)
        
    if label:
        labels = {"o": [0,0,0],
                  "a": vecs[0,:],
                  "b": vecs[1,:],
                  "c": vecs[2,:]}
        for key in labels.keys():
            text, pos = key, labels[key]
            if viewer is None:
                view.addLabel(text, {"position": {"x": pos[0], "y": pos[1], "z": pos[2]}, 
                                     "fontColor":"black",
                                     "backgroundColor": "white", 
                                     "fontsize": 12,
                                     "backgroundOpacity": "0.1",
                                    })
            else:
                view.addLabel(text, {"position": {"x": pos[0], "y": pos[1], "z": pos[2]}, 
                                     "fontColor":"black",
                                     "backgroundColor": "white", 
                                     "fontsize": 12,
                                     "backgroundOpacity": "0.1",
                                    },
                              viewer = viewer)
        
def addlines(view, orig, axes, viewer=None):
    for ax in axes:
        ax += orig
        if viewer is None:
            view.addArrow({"start": {"x":orig[0], "y": orig[1], "z":orig[2]}, 
                       "end": {"x":ax[0], "y":ax[1], "z":ax[2]}})
        else:
            view.addArrow({"start": {"x":orig[0], "y": orig[1], "z":orig[2]}, 
                           "end": {"x":ax[0], "y":ax[1], "z":ax[2]}},
                           viewer=viewer)



def display_atomic_crystal(structure, height=300, width=600, scale=0.25, radius=0.10, show_wp=True):
    view = py3Dmol.view(height=height, width=width)
    if structure.dim == 0:
        fmt = "xyz"
    else:
        fmt = "cif"

    txt = structure.to_file(fmt=fmt)
    view.addModel(txt, fmt, {'doAssembly':True,'duplicateAssemblyAtoms':True})
    view.setStyle({'sphere':{'colorscheme':'Jmol','scale':scale},
                   'stick':{'colorscheme':'Jmol', 'radius':radius}})
    if structure.dim != 0:
        view.addUnitCell()
        if show_wp:
            view.setStyle({'sym':2},{'sphere':{'scale':scale*1.1,'color':'blue'}})
    return view.zoomTo()

def display_molecular_crystal(structure, height=300, width=600, axis=False, axis_label=True):
    view = py3Dmol.view(height=height, width=width)
    for site in structure.mol_sites:
        for i in range(site.multiplicity):
            mol = site.get_mol_object(i)
            view.addModel(mol.to(fmt='xyz'), 'xyz')
            if axis:
                axes = site.get_principle_axes(mol.cart_coords)
                addlines(view, np.mean(mol.cart_coords, axis=0), axes.T*5)

    addBox(view, structure.lattice_matrix, axis_label)
    view.setStyle({'stick':{'colorscheme':'greenCarbon'}})

    return view.zoomTo()

def display_molecular_site(site, id=None, height=300, width=800, axis=True, axis_label=True):
    view = py3Dmol.view(height=height, width=width, linked=False, viewergrid=(1,2))

    view.addModel(site.mol.to(fmt='xyz'), 'xyz', viewer=(0,0))
    view.setStyle({'stick':{'colorscheme':'blueCarbon'}}, viewer=(0,0))

    if id is not None:
        ids = range(site.multiplicity)
    else:
        ids = [id]

    for i in range(site.multiplicity):
        mol = site.get_mol_object(i)
        view.addModel(mol.to(fmt='xyz'), 'xyz', viewer=(0,1))
        if axis:
            axes = site.get_principle_axes(mol.cart_coords)
            addlines(view, np.mean(mol.cart_coords, axis=0), axes.T*5, viewer=(0,1))
    
    addBox(view, site.lattice, axis_label, viewer=(0,1))
    view.setStyle({'stick':{'colorscheme':'greenCarbon'}})

    return view.zoomTo()


