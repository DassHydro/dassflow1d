import argparse
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dassflow1d import write_mesh
import dassflow1d.m_mesh as m_mesh

def load_segments_shapefile(segments_file):

    # Load segments shapefile
    print("Load segments shapefile")
    segments = gpd.read_file(segments_file)
    print("- Number of segments: %i" % segments.index.size)

    # Find ID attribute
    if "reach_id" in segments.columns:
        # Use 'reach_id' attribute (SWORD)
        id_attribute = "reach_id"
    elif "seg_id" in segments.columns:
        # Use 'seg_id' attribute (Effective models)
        id_attribute = "seg_id"
    elif "id" in segments.columns:
        # Use default attribute 'id'
        id_attribute = "id"
    else:
        raise RuntimeError("No attribute found for segment ID (possible attributes are: reach_id, seg_id, id)")
    print("- Using ID attribute: %s" % id_attribute)

    # Find upstream connectivity attribute
    if "rch_id_up" in segments.columns:
        # Use 'rch_id_up' attribute (SWORD)
        us_conn_attribute = "rch_id_up"
    elif "up_segs" in segments.columns:
        # Use 'up_segs' attribute (Effective models)
        us_conn_attribute = "up_segs"
    elif "us_conn" in segments.columns:
        # Use default attribute 'us_conn'
        us_conn_attribute = "us_conn"
    else:
        raise RuntimeError("No attribute found for upstream connectivity (possible attributes are: rch_id_up, up_segs, us_conn)")
    print("- Using upstream connectivity attribute: %s" % us_conn_attribute)

    # Find downstream connectivity attribute
    if "rch_id_dn" in segments.columns:
        # Use 'up_segs' attribute (SWORD)
        ds_conn_attribute = "rch_id_dn"
    elif "down_seg" in segments.columns:
        # Use 'down_seg' attribute (Effective models, dentritic networks)
        ds_conn_attribute = "down_seg"
    elif "down_segs" in segments.columns:
        # Use 'down_segs' attribute (Effective models)
        ds_conn_attribute = "down_segs"
    elif "ds_conn" in segments.columns:
        # Use default attribute 'ds_conn'
        ds_conn_attribute = "ds_conn"
    else:
        raise RuntimeError("No attribute found for downstream connectivity (possible attributes are: rch_id_dn, down_seg, down_segs, ds_conn)")
    print("- Using downstream connectivity attribute: %s" % ds_conn_attribute)

    return {"data": segments.loc[:, [id_attribute, us_conn_attribute, ds_conn_attribute, "geometry"]],
            "id_attribute": id_attribute,
            "us_conn_attribute": us_conn_attribute,
            "ds_conn_attribute": ds_conn_attribute}


def load_sections_shapefile(sections_file):

    # Load sections shapefile
    print("Load sections shapefile")
    sections = gpd.read_file(sections_file)
    print("- Number of sections: %i" % sections.index.size)

    # Find ID attribute
    if "node_id" in sections.columns:
        # Use 'node_id' attribute (SWORD)
        id_attribute = "node_id"
    elif "vs_id" in sections.columns:
        # Use 'vs_xs' attribute (Effective models)
        id_attribute = "vs_id"
    elif "id" in sections.columns:
        # Use default attribute 'id'
        id_attribute = "id"
    else:
        raise RuntimeError("No attribute found for section ID (possible attributes are: node_id, vs_id, id)")
    print("- Using ID attribute: %s" % id_attribute)

    # Find segments ID attribute
    if "reach_id" in sections.columns:
        # Use 'reach_id' attribute (SWORD)
        segment_attribute = "reach_id"
    elif "segment" in sections.columns:
        # Use 'segment' attribute (Effective models)
        segment_attribute = "segment"
    elif "seg_id" in sections.columns:
        # Use default attribute 'seg_id'
        segment_attribute = "seg_id"
    else:
        raise RuntimeError("No attribute found for segment (possible attributes are: reach_id, segment, seg_id)")
    print("- Using segment attribute: %s" % segment_attribute)

    # Find curvilinear abscissa attribute
    if "dist_out" in sections.columns:
        # Use 'dist_out' attribute (SWORD)
        xs_attribute = "dist_out"
    elif "vs_xs" in sections.columns:
        # Use 'vs_xs' attribute (Effective models)
        xs_attribute = "vs_xs"
    elif "xs" in sections.columns:
        # Use default attribute 'xs'
        xs_attribute = "xs"
    else:
        xs_attribute = None
    if xs_attribute is None:
        print("- No curvilinear abscissa attribute found, curvilinear abscissa will be computed from geometries")
    else:
        print("- Using curvilinear abscissa attribute: %s" % xs_attribute)

    # Find bathymetry attribute
    if "vs_bathy" in sections.columns:
        # Use 'vs_bathy' attribute (Effective models)
        bathy_attribute = "vs_bathy"
    elif "bathy" in sections.columns:
        # Use default attribute 'bathy'
        bathy_attribute = "bathy"
    else:
        bathy_attribute = None
    if bathy_attribute is None:
        print("- No bathy attribute found, bathy will be estimated")
    else:
        print("- Using bathy attribute: %s" % xs_attribute)

    sections[id_attribute] = np.round(sections[id_attribute]).astype(int)

    return {"data": sections.loc[:, [id_attribute, segment_attribute, xs_attribute, "geometry"]],
            "id_attribute": id_attribute,
            "segment_attribute": segment_attribute,
            "xs_attribute": xs_attribute,
            "bathy_attribute": bathy_attribute}


def load_profiles_csv(profiles_file, section_attribute, sep=";"):

    print("Load sections profiles")
    profiles = pd.read_csv(profiles_file, sep=sep)
    print(profiles.columns)

    # Find H attribute
    if "wse" in profiles.columns:
        # Use 'wse' attribute (SWOT)
        H_attribute = "wse"
    elif "elev" in profiles.columns:
        # Use 'elev' attribute (Effective models)
        H_attribute = "elev"
    elif "H" in profiles.columns:
        # Use default attribute 'H'
        H_attribute = "H"
    else:
        raise RuntimeError("No attribute found for H (possible attributes are: wse, elev, H)")
    print("- Using H attribute: %s" % H_attribute)

    # Find W attribute
    if "width" in profiles.columns:
        # Use 'width' attribute (SWOT)
        W_attribute = "width"
    elif "W" in profiles.columns:
        # Use default attribute 'W'
        W_attribute = "W"
    else:
        raise RuntimeError("No attribute found for W (possible attributes are: width, W)")
    print("- Using W attribute: %s" % W_attribute)

    profiles[section_attribute] = np.round(profiles[section_attribute]).astype(int)
    print(profiles[section_attribute])

    # Check section attribute is in columns
    if not section_attribute in profiles.columns:
        raise RuntimeError("Section attribute '%s' not found in columns of profiles table")

    return {"data": profiles.loc[:, [section_attribute, H_attribute, W_attribute]],
            "H_attribute": H_attribute,
            "W_attribute": W_attribute}


def create_dassflow1d_mesh(segments_file, sections_file, profiles_file, mesh_file, profiles_sep=";"):

    # Load segments
    segments_data = load_segments_shapefile(segments_file)

    # Load sections
    sections_data = load_sections_shapefile(sections_file)

    # Load profile CSV file
    profiles_data = load_profiles_csv(profiles_file, sections_data["id_attribute"], profiles_sep)

    # Initialize mesh
    print("Create mesh")
    segments_df = segments_data["data"]
    sections_df = sections_data["data"]
    print(sections_df["segment"].unique())
    choice = input()
    profiles_df = profiles_data["data"]
    mesh = m_mesh.Mesh(sections_df.index.size, segments_df.index.size)
    print("- Initialise mesh with %i sections and %i segments" % (sections_df.index.size, segments_df.index.size))

    # Setup sections
    sections_offset = 0
    print("- Setup sections")
    for iseg, seg_index in enumerate(segments_df.index):
        # Retrieve segments sections
        seg_id = segments_df.loc[seg_index, segments_data["id_attribute"]]
        segments_sections = sections_df[sections_df[sections_data["segment_attribute"]] == seg_id]
        print("- Setup %i (%i-%i) sections for segment %i (ID=%i)" % (segments_sections.index.size, sections_offset, sections_offset+segments_sections.index.size, seg_index, iseg))

        if sections_data["xs_attribute"] is None:
            raise NotImplementedError("Computing curvilinear abscissae from geometry is not implemented yet")

        # Sort segments sections by decreasing curvilinear abscissa
        segments_sections = segments_sections.sort_values(by=sections_data["xs_attribute"], ascending=False)

        # Setup sections using profile table
        for ics, index in enumerate(segments_sections.index):
            section_profiles_df = profiles_df[profiles_df[sections_data["id_attribute"]] == segments_sections.loc[index, sections_data["id_attribute"]]]
            section_df = sections_df[sections_df[sections_data["id_attribute"]] == segments_sections.loc[index, sections_data["id_attribute"]]]
            xs = section_df.loc[:, sections_data["xs_attribute"]].iloc[0]
            H = section_profiles_df.loc[:, profiles_data["H_attribute"]].values
            W = section_profiles_df.loc[:, profiles_data["W_attribute"]].values
            print(section_df.loc[:, sections_data["id_attribute"]].iloc[0], ics, H, W)
            isort = np.argsort(H)
            Hs = H[isort]
            Ws = W[isort]
            if sections_data["bathy_attribute"] is None:
                bathy = Hs[0] - 5.0
            else:
                bathy = section_df.loc[:, sections_data["bathy_attribute"]].iloc[0]
            # TODO Check that Ws is monotonicaly increasing
            # mesh.setup_crosssection(sections_offset+ics+1, Hs.size)
            mesh.cs[sections_offset+ics].x = xs
            mesh.cs[sections_offset+ics].set_coords(xs, 0.0)
            mesh.cs[sections_offset+ics].set_levels(Hs, Ws)
            mesh.cs[sections_offset+ics].bathy = bathy
        sections_offset += segments_sections.index.size
            # if np.any(OBc == 1):
            #     OB_indices = np.ravel(np.argwhere(OBc == 1))
            #     mesh.cs[sections_offset+ics].set_overbanks_levels(OB_indices[0], OB_indices[0])            

    # TODO Ensure segments are ordered from upstream to downstream

    # Setup segments
    segments_df = segments_data["data"]
    sections_df = sections_data["data"]
    sections_offset = 0
    for iseg, index in enumerate(segments_df.index):

        seg_id = segments_df.loc[index, segments_data["id_attribute"]]

        # Compute list of upstream segments indices
        us_conn = segments_df.loc[index, segments_data["us_conn_attribute"]]
        if us_conn is None:
            us_indices = [-1]
        else:
            us_ids = [int(id_str) for id_str in us_conn.split(",")]
            us_indices = []
            for us_id in us_ids:
                us_row = segments_df[segments_df[segments_data["id_attribute"]] == us_id]
                if us_row.index.size == 0:
                    raise RuntimeError("Upstream segment with ID %i not found in segments shapefile" % us_id)
                elif us_row.index.size > 1:
                    raise RuntimeError("Upstream segment with ID %i found multiple times in segments shapefile" % us_id)
                us_indices.append(us_row.index[0])
        print("us_conn[%i]:" % iseg, us_conn, type(us_conn))

        # Compute list of downstream segments indices
        ds_conn = segments_df.loc[index, segments_data["ds_conn_attribute"]]
        if ds_conn is None:
            ds_indices = [-1]
        else:
            if isinstance(ds_conn, str):
                ds_ids = [int(id_str) for id_str in ds_conn.split(",")]
            else:
                ds_ids = [ds_conn]
            ds_indices = []
            for ds_id in ds_ids:
                if ds_id == -9 or ds_id == -1:
                    ds_indices.append(-1)
                else:
                    ds_row = segments_df[segments_df[segments_data["id_attribute"]] == ds_id]
                    if ds_row.index.size == 0:
                        raise RuntimeError("Downstream segment with ID %i not found in segments shapefile" % ds_id)
                    elif ds_row.index.size > 1:
                        raise RuntimeError("Downstream segment with ID %i found multiple times in segments shapefile" % ds_id)
                    ds_indices.append(ds_row.index[0])
        print("ds_conn[%i]:" % iseg, ds_conn, type(ds_conn))

        # Compute list of sections
        segments_sections = sections_df[sections_df[sections_data["segment_attribute"]] == seg_id]
        mesh.setup_segment(iseg, sections_offset, segments_sections.index.size, us_indices, ds_indices[0])

    mesh.check()

    write_mesh(mesh_file, mesh)

    # with open(geometry_file, "r") as geometry_fp:
    #     content = geometry_fp.readlines()

    # line_index = 0
    # profile_y = []
    # profile_z = []
    # profile_flag = []
    # segments = []
    # current_segment_name = None
    # sections_count = 0
    # for line_index, line in enumerate(content):

    #     words = line.strip().split(" ")

    #     if words[0] == "PROFIL":

    #         # Store previous profile
    #         if len(profile_y) > 0:
    #             segments[-1]["profiles"].append({"xs": xs,
    #                                              "name": section_name,
    #                                              "y": np.array(profile_y),
    #                                              "z": np.array(profile_z),
    #                                              "flag": np.array(profile_flag)})
    #             sections_count += 1

    #         # Start reading profile
    #         xs = float(words[3])
    #         segment_name = words[1]
    #         section_name = words[2]
    #         if segment_name != current_segment_name:
    #             segments.append({"name": segment_name,
    #                              "profiles": []})
    #         current_segment_name = segment_name
    #         profile_y = []
    #         profile_z = []
    #         profile_flag = []

    #     else:

    #         # Read profile data
    #         profile_y.append(float(words[0]))
    #         profile_z.append(float(words[1]))
    #         profile_flag.append(words[2])

    # # Store last profile
    # if len(profile_y) > 0:
    #     segments[-1]["profiles"].append({"xs": xs,
    #                                      "name": section_name,
    #                                      "y": np.array(profile_y),
    #                                      "z": np.array(profile_z),
    #                                      "flag": np.array(profile_flag)})
    # sections_count += 1
        
    # if len(segments) > 1:
    #     raise NotImplementedError("Converting geometry with multiples segments is not implemented yet.")

    # # Initialise mesh
    # mesh = m_mesh.Mesh(sections_count, len(segments))

    # sections_offset = 0
    # for iseg, segment in enumerate(segments):

    #     # Setup segment
    #     mesh.setup_segment(iseg, sections_offset, len(segment["profiles"]), -1, -1)

    #     # Compute trapezoidal section
    #     for ics, profile in enumerate(segment["profiles"]):
    #         # plt.plot(profile["y"], profile["z"])
    #         # plt.plot(profile["y"][profile["flag"] == "T"], profile["z"][profile["flag"] == "T"], 'r.')
    #         # plt.title("Segment %s, profile %s (xs=%f)" % (segment["name"], profile["name"], profile["xs"]))
    #         # plt.show()

    #         y = profile["y"]
    #         z = profile["z"]
    #         flag = profile["flag"]
    #         indices_sorted = np.argsort(z)
    #         H = np.zeros(z.size)
    #         W = np.zeros(z.size)
    #         OB = np.zeros(z.size, dtype=int)
    #         for i in range(0, len(indices_sorted)):
    #             ww = 0
    #             zw = z[indices_sorted[i]]
    #             # print(i, zw)
    #             for j in range(0, len(indices_sorted)-1):
    #                 if z[j] <= zw and z[j+1] <= zw:
    #                     ww += y[j+1] - y[j]
    #             H[i] = zw
    #             W[i] = ww
    #             OB[i] = flag[indices_sorted[i]] == "T"

    #         # plt.plot(W, H)
    #         # plt.plot(W[OB == 1], H[OB == 1], "r.")
    #         # plt.show()

    #         # Remove doublons
    #         i = 2
    #         Hc = []
    #         Wc = []
    #         OBc = []
    #         if sections_offset+ics == 73:
    #             print("Bottom0 %i %f %f" % (0, H[0], W[0]))
    #             print("Bottom1 %i %f %f" % (1, H[1], W[1]))
    #         while i < len(H):
    #             j = i
    #             while H[i] - H[j-1] < 0.00001 or W[i] - W[j-1] < 0.00001:
    #                 if sections_offset+ics == 73:
    #                     print("Doublon %i %f %f" % (i, H[i], W[i]))
    #                 if i == len(H)-1:
    #                     break
    #                 i = i + 1
    #             if i == len(H)-1:
    #                 if i == j and H[i] - H[j-1] > 0.00001 and W[i] - W[j-1] > 0.00001:
    #                     if sections_offset+ics == 73:
    #                         print("Append* %i %f %f" % (i, H[i], W[i]))
    #                     Hc.append(H[i])
    #                     Wc.append(W[i])
    #                     OBc.append(OB[i])
    #                 break
    #             else:
    #                 if sections_offset+ics == 73:
    #                     print("Append+ %i (%i) %f %f (%i %f %f)" % (i, j, H[i], W[i], i-1, H[i-1], W[i-1]))
    #                 Hc.append(H[i])
    #                 Wc.append(W[i])
    #                 OBc.append(OB[i])
    #                 i = i + 1
    #             # if i == j:
    #             #     i = i + 1

    #         A0 = 0.5 * (Hc[1] - H[0]) * Wc[1]
    #         bathy = Hc[1] - A0 / Wc[1]

    #         # Setup section
    #         print("ICS:", sections_offset+ics, len(Hc))
    #         # mesh.setup_crosssection(sections_offset+ics, H.size-1)
    #         mesh.cs[sections_offset+ics].x = profile["xs"]
    #         mesh.cs[sections_offset+ics].set_coords(profile["xs"], 0.0)
    #         mesh.cs[sections_offset+ics].set_levels(Hc[1:], Wc[1:])
    #         mesh.cs[sections_offset+ics].bathy = bathy
    #         if np.any(OBc == 1):
    #             OB_indices = np.ravel(np.argwhere(OBc == 1))
    #             mesh.cs[sections_offset+ics].set_overbanks_levels(OB_indices[0], OB_indices[0])



    # write_mesh("mesh.geo", mesh)



if __name__ == "__main__":

    parser = argparse.ArgumentParser("Create DassFlow-1D mesh from definition files")
    parser.add_argument("segments_file", type=str, help="Path to the shapefile with segments definitions")
    parser.add_argument("sections_file", type=str, help="Path to the shapefile with sections definitions")
    parser.add_argument("profiles_files", type=str, help="Path to the CSV file with the table defining (y,z) profiles")
    parser.add_argument("-o", dest="mesh_file", type=str, default="mesh.geo", help="Path to the output mesh file")
    parser.add_argument("-sep", type=str, default=";", help="Separator for the CSV file")
    #parser.add_argument("geometry_file", type=str, help="Geometry file")

    args = parser.parse_args()

    # Convert geometry file to DassFlow-1D mesh
    create_dassflow1d_mesh(args.segments_file, args.sections_file, args.profiles_files, args.mesh_file, profiles_sep=args.sep)
