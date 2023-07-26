import numpy as np
import matplotlib.pyplot as plt


def restrict_between_levees(t, z, ibanks, levee_threshold=1.0):
    """ Restrict (t,z) values between levees
    
        Parameters:
            t (numpy.ndarray): array of transversal distance from left to right
            z (numpy.ndarray): array of elevations corresponding to values in t
            ibanks (tuple): indices (in array t) of banks
            levee_threshold (float): Threshold for removing the points behind levees
            fname (str): Path to the plot file
        Return:
            t (numpy.ndarray): array of transversal distance restricted tosubset between levees
            z (numpy.ndarray): array of elevations restricted tosubset between levees
            ibanks (numpy.ndarray): Updated indices (in array t) of banks
    """

    # Compute restriction on left side
    dz = z[1:ibanks[0]+1] - z[0:ibanks[0]]
    indices = np.where(dz > levee_threshold)[0]
    if len(indices) == 0:
        ilevee_left = 0
    else:
        ilevee_left = indices[-1]+1

    # Compute restriction on right side
    dz = z[ibanks[1]:-1] - z[ibanks[1]+1:]
    indices = np.where(dz > levee_threshold)[0]
    if len(indices) == 0:
        ilevee_right = t.size-1
    else:
        ilevee_right = indices[0] + ibanks[1]
    
    # Compute new arrays t and z
    t = t[ilevee_left:ilevee_right+1]
    z = z[ilevee_left:ilevee_right+1]

    # Update banks indices
    if ilevee_left > 0:
        ibanks[0] -= ilevee_left
        ibanks[1] -= ilevee_left
        
    return t, z, ibanks


def effective_section(t, z, ibanks, hmin=0.1, levee_threshold=1.0, fname=None):
    """ Compute a effective section from a real section (t,z(t)).
    
        Parameters:
            t (numpy.ndarray): array of transversal distance from left to right
            z (numpy.ndarray): array of elevations corresponding to values in t
            ibanks (tuple): Indices (in array t) of banks
            hmin (float): Minimal depth (i.e. unobseved depth)
            levee_threshold (float): Threshold for removing the points behind levees
            fname (str): Path to the plot file
        Return:
            bathy (float): bathymetry elevation
            H (numpy.ndarray): array of elevations
            A (numpy.ndarray): array of flow areas
            W (numpy.ndarray): array of top widths
    """
    
    # Initialise figure
    if fname is not None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
        ax1.plot([t[ibanks[0]], t[ibanks[1]]], [z[ibanks[0]], z[ibanks[1]]], "g+")
        ax1.plot(t, z, "k-")
        
    # Restrict between levees
    if levee_threshold is not None:
        t, z, ibanks = restrict_between_levees(t, z, ibanks, levee_threshold)
    
    # Show banks/levees points
    if fname is not None:
        ax1.plot(t, z, "b-")
        ax1.plot([t[ibanks[0]], t[ibanks[1]]], [z[ibanks[0]], z[ibanks[1]]], "ro")
    
    # Sort elevations
    zp = np.sort(z)
    
    # Compute reference elevation (to compute areas using integrals)
    z0 = np.minimum(0.0, zp[0] - 1.0)
    
    # Initialise lists of H, W and A
    H = []
    W = []
    A = []
    
    # Compute highest and lowest elevations between adjacent tuple of (t,z)
    zup = np.maximum(z[0:-1], z[1:])
    zdn = np.minimum(z[0:-1], z[1:])
    
    # Compute dW (width between to values of t)
    dW = t[1:] - t[0:-1]

    # Initialise dA (flow area between to values of t)
    dA = np.zeros(zup.size)
    for i in range(0, zp.size):
        
        # Compute dA = area between bathymetry and zp[i]
        dA = np.maximum(0.0, (zp[i] - z0) * dW - (0.5*(zup + zdn) - z0) * dW)
        
        # Compute Ap=sum(dA)
        Ap = np.sum(dA)

        # Compute Wp=sum(dW where dA > 0.0)
        Wp = np.sum((dA > 0.0) * dW)

        if zp[i] - zp[0] >= hmin and np.any(dA > 0.0):
            
            if len(A) == 0:
                
                # Compute W0 = Wp
                W0 = Wp
                
                # Compute bathy
                bathy = zp[i] - Ap / W0
                
                # Initial set of H, A and W
                H.append(zp[i])
                A.append(Ap)
                W.append(Wp)

            elif Ap >= A[-1] + 0.1:
                
                # Append new set of H, A and W
                H.append(zp[i])
                A.append(Ap)
                W.append(Wp)
                
    # Update and save figure
    if fname is not None:
        ax2.plot(W, H)
        plt.savefig(fname)
        plt.close(fig)
        
    return bathy, H, A, W

def save_mesh(mesh, fname, version=1.0, shape_model="linear", proj="WGS84"):

    # Open file
    fout = open(fname, "w")
    
    # Write header
    fout.write("DASSFLOW-1D MESH V%.2i\n" % version)
    fout.write("# Header\n")
    fout.write("%i %i %s" % (mesh.ncs, mesh.nseg, shape_model))
    if version >= 1.3:
        fout.write("%s\n" % proj)
    else:
        fout.write("\n")
        
    # Write sections
    fout.write("# Sections\n")
    for ics in range(0, mesh.ncs):
        
        # Write global parameters and number of levels
        fout.write("%i" % (ics+1))
        if version >= 1.3:
            fout.write(" %f" % mesh.cs[ics].x)
        fout.write(" %f %f %f %i\n" % (mesh.cs[ics].coord.x, mesh.cs[ics].coord.y, mesh.cs[ics].bathy, 
                                            mesh.cs[ics].nlevels))
        
        # Write levels
        for ilevel in range(0, mesh.cs[ics].nlevels):
            fout.write("%f %f 0.0" % (mesh.cs[ics].level_heights[ilevel], mesh.cs[ics].level_widths[ilevel]))
            if version >= 1.0:
                if ilevel+1 >= mesh.cs[ics].ob_levels[0]:
                    flag = 10
                else:
                    flag = 0
                if ilevel+1 >= mesh.cs[ics].ob_levels[1]:
                    flag += 1
                fout.write(" %02i\n" % flag)
            else:
                fout.write("\n")
                
    # Write segments
    fout.write("# Segments\n")
    for iseg in range(0, mesh.nseg):
        fout.write("%i %i %i %i\n" % (iseg+1, mesh.seg[iseg].first_cs, mesh.seg[iseg].last_cs, 
                                      len(mesh.seg[iseg].us_seg)))
        for i in range(0, len(mesh.seg[iseg].us_seg)):
            fout.write(" %i" % mesh.seg[iseg].us_seg[i])
        fout.write(" %i\n" % mesh.seg[iseg].ds_seg)
        
    # Close file
    fout.close()
        
        
                
    
