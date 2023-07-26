# coding: utf-8
from __future__ import print_function

import argparse
import numpy as np
from dassflow1d.utils.mesh_utils import effective_section, save_mesh
from dassflow1d.utils.hecras_geoparser import HecRasGeoParser
from dassflow1d.utils.stdout_utils import stdout_progress
from dassflow1d import m_mesh


def convert_hecras_geofile_to_mesh(geofile, meshfile, hmin=0.1, levee_threshold=None, plot_prefix=None, 
                                   real_sections=True, debug=False):
    """ Restrict (t,z) values between levees
    
        Parameters:
            geofile (str): Path to the HEC-RAS geo file
            meshfile (str): Path to the output mesh file
            hmin (float): Minimal depth (i.e. unobseved depth)
            levee_threshold (float): Threshold for removing the points behind levees
            plot_prefix (str): Prefix for the plot files
            real_sections (bool): True if only real (not interpolated) sections are converted
            debug (bool): True for debug mode
    """
    
    # TODO implement debug mode
    # TODO implement multi-segments
  
    #------------------------------------------------------------------------------------------------------------------
    # Load HEC RAS geometry
    #------------------------------------------------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print(" LOADING HEC-RAS GEOMETRY")
    print("=" * 80 + "\n")
    parser = HecRasGeoParser(geofile, parse_only_real_sections=real_sections)
    print("Number of sections: %i" % len(parser.sections))
    print("Number of reaches: %i" % len(parser.reaches))

    #------------------------------------------------------------------------------------------------------------------
    # Retreave single reach
    #------------------------------------------------------------------------------------------------------------------
    if len(parser.reaches) > 1:
        raise NotImplementedError("Converting geo with multiple reaches is not implemented yet")
    reach = parser.reaches[0]
    
    #------------------------------------------------------------------------------------------------------------------
    # Convert to mesh
    #------------------------------------------------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print(" CREATING MESH")
    print("=" * 80 + "\n")
    # Initialise mesh
    mesh = m_mesh.Mesh(len(reach.sections), 1)
    
    # Setup cross-sections
    for ics in range(0, len(reach.sections)):
        
        stdout_progress("Progress:", ics, maxval=len(reach.sections))
        # Retrieve t and z arrays and indices of banks
        t = reach.sections[ics].t
        z = reach.sections[ics].z
        ibanks = reach.sections[ics].ibanks

        # Compute effective section
        if plot_prefix is not None:
            plot_file = "%s_reach00_xs%04i.png" % (plot_prefix, ics+1) 
        else:
            plot_file = None
        bathy, H, A, W = effective_section(t, z, ibanks, hmin, levee_threshold, plot_file)
        
        mesh.cs[ics].set_levels(H, W)
        ibank_left = np.argmin(H >= z[ibanks[0]])
        ibank_right = np.argmin(H >= z[ibanks[1]])
        mesh.cs[ics].set_overbanks_levels(ibank_left, ibank_right)
        mesh.cs[ics].bathy = bathy
        mesh.cs[ics].coord.x = 0.0
        mesh.cs[ics].coord.y = 0.0

    stdout_progress("Progress:", len(reach.sections), maxval=len(reach.sections))
    
    # Setup segments
    mesh.seg[0].set_crosssections_range(0, len(reach.sections))
    mesh.seg[0].set_connectivity([-1], -1)
    
    save_mesh(mesh, meshfile, version=1.3, proj="None")

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Generate priors using Metroman method")
    parser.add_argument("geofile", type=str,
                        help="Path to the HEC-RAS geo file")
    parser.add_argument("meshfile", type=str,
                        help="Path to the output mesh file")
    parser.add_argument("-hmin", type=float, default=0.1, 
                        help="Minimal (unobserved) depth")
    parser.add_argument("-levee_threshold", type=float, default=None,
                        help="Threshold for restriction between levees")
    parser.add_argument("-plot-prefix", type=str, default=None, 
                        help="Prefix for plot files")
    parser.add_argument("--real-sections", dest="debug_mode", action='store_true',
                        help="Converts only real (not interpolated) sections")
    
    args = parser.parse_args()
    convert_hecras_geofile_to_mesh(args.geofile, args.meshfile, args.hmin, args.levee_threshold, args.plot_prefix)
    
