import argparse
import matplotlib.pyplot as plt
import numpy as np

from dassflow1d import write_mesh
import dassflow1d.m_mesh as m_mesh

def create_dassflow1d_mesh(geometry_file):

    with open(geometry_file, "r") as geometry_fp:
        content = geometry_fp.readlines()

    line_index = 0
    profile_y = []
    profile_z = []
    profile_flag = []
    segments = []
    current_segment_name = None
    sections_count = 0
    for line_index, line in enumerate(content):

        words = line.strip().split(" ")

        if words[0] == "PROFIL":

            # Store previous profile
            if len(profile_y) > 0:
                segments[-1]["profiles"].append({"xs": xs,
                                                 "name": section_name,
                                                 "y": np.array(profile_y),
                                                 "z": np.array(profile_z),
                                                 "flag": np.array(profile_flag)})
                sections_count += 1

            # Start reading profile
            xs = float(words[3])
            segment_name = words[1]
            section_name = words[2]
            if segment_name != current_segment_name:
                segments.append({"name": segment_name,
                                 "profiles": []})
            current_segment_name = segment_name
            profile_y = []
            profile_z = []
            profile_flag = []

        else:

            # Read profile data
            profile_y.append(float(words[0]))
            profile_z.append(float(words[1]))
            profile_flag.append(words[2])

    # Store last profile
    if len(profile_y) > 0:
        segments[-1]["profiles"].append({"xs": xs,
                                         "name": section_name,
                                         "y": np.array(profile_y),
                                         "z": np.array(profile_z),
                                         "flag": np.array(profile_flag)})
    sections_count += 1
        
    if len(segments) > 1:
        raise NotImplementedError("Converting geometry with multiples segments is not implemented yet.")

    # Initialise mesh
    mesh = m_mesh.Mesh(sections_count, len(segments))

    sections_offset = 0
    for iseg, segment in enumerate(segments):

        # Setup segment
        mesh.setup_segment(iseg, sections_offset, len(segment["profiles"]), -1, -1)

        # Compute trapezoidal section
        for ics, profile in enumerate(segment["profiles"]):
            # plt.plot(profile["y"], profile["z"])
            # plt.plot(profile["y"][profile["flag"] == "T"], profile["z"][profile["flag"] == "T"], 'r.')
            # plt.title("Segment %s, profile %s (xs=%f)" % (segment["name"], profile["name"], profile["xs"]))
            # plt.show()

            y = profile["y"]
            z = profile["z"]
            flag = profile["flag"]
            indices_sorted = np.argsort(z)
            H = np.zeros(z.size)
            W = np.zeros(z.size)
            OB = np.zeros(z.size, dtype=int)
            for i in range(0, len(indices_sorted)):
                ww = 0
                zw = z[indices_sorted[i]]
                # print(i, zw)
                for j in range(0, len(indices_sorted)-1):
                    if z[j] <= zw and z[j+1] <= zw:
                        ww += y[j+1] - y[j]
                H[i] = zw
                W[i] = ww
                OB[i] = flag[indices_sorted[i]] == "T"

            # plt.plot(W, H)
            # plt.plot(W[OB == 1], H[OB == 1], "r.")
            # plt.show()

            # Remove doublons
            i = 2
            Hc = []
            Wc = []
            OBc = []
            if sections_offset+ics == 73:
                print("Bottom0 %i %f %f" % (0, H[0], W[0]))
                print("Bottom1 %i %f %f" % (1, H[1], W[1]))
            while i < len(H):
                j = i
                while H[i] - H[j-1] < 0.00001 or W[i] - W[j-1] < 0.00001:
                    if sections_offset+ics == 73:
                        print("Doublon %i %f %f" % (i, H[i], W[i]))
                    if i == len(H)-1:
                        break
                    i = i + 1
                if i == len(H)-1:
                    if i == j and H[i] - H[j-1] > 0.00001 and W[i] - W[j-1] > 0.00001:
                        if sections_offset+ics == 73:
                            print("Append* %i %f %f" % (i, H[i], W[i]))
                        Hc.append(H[i])
                        Wc.append(W[i])
                        OBc.append(OB[i])
                    break
                else:
                    if sections_offset+ics == 73:
                        print("Append+ %i (%i) %f %f (%i %f %f)" % (i, j, H[i], W[i], i-1, H[i-1], W[i-1]))
                    Hc.append(H[i])
                    Wc.append(W[i])
                    OBc.append(OB[i])
                    i = i + 1
                # if i == j:
                #     i = i + 1

            A0 = 0.5 * (Hc[1] - H[0]) * Wc[1]
            bathy = Hc[1] - A0 / Wc[1]

            # Setup section
            print("ICS:", sections_offset+ics, len(Hc))
            # mesh.setup_crosssection(sections_offset+ics, H.size-1)
            mesh.cs[sections_offset+ics].x = profile["xs"]
            mesh.cs[sections_offset+ics].set_coords(profile["xs"], 0.0)
            mesh.cs[sections_offset+ics].set_levels(Hc[1:], Wc[1:])
            mesh.cs[sections_offset+ics].bathy = bathy
            if np.any(OBc == 1):
                OB_indices = np.ravel(np.argwhere(OBc == 1))
                mesh.cs[sections_offset+ics].set_overbanks_levels(OB_indices[0], OB_indices[0])



    write_mesh("mesh.geo", mesh)






if __name__ == "__main__":

    parser = argparse.ArgumentParser("Create DassFlow-1D mesh from definition files")
    parser.add_argument("segments_file", type=str, help="Path to the shapefile with segments definitions")
    parser.add_argument("sections_file", type=str, help="Path to the shapefile with sections definitions")
    parser.add_argument("profiles_files", type=str, help="Path to the CSV file with the table defining (y,z) profiles")
    #parser.add_argument("geometry_file", type=str, help="Geometry file")

    args = parser.parse_args()

    # Convert geometry file to DassFlow-1D mesh
    create_dassflow1d_mesh(args.segments_file, args.sections_file, args.profiles_files)
