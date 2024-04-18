import argparse
import glob
import os


def change_version_in_file(fname, version, comment_character):

    rootdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Read source file
    with open(os.path.join(rootdir, fname), "r") as src:
        content = src.readlines()

    # Check third line content version
    words = content[2].split("Version")
    expected = "%s                    DassFlow1D " % comment_character
    expected_variant = "%s                    DassFlow " % comment_character
    if not words[0] == expected and not words[0] == expected_variant:
        print("words[0]='%s'" % words[0])
        raise ValueError("Third line in file %s must be %sVersion X.X" % (fname, expected))
    content[2] = "%sVersion %s\n" % (expected, version)

    # Modify file
    with open(os.path.join(rootdir, fname), "w") as dst:
        for line in content:
            dst.write(line)

    print("Modified file %s" % fname)


def change_version_in_makefiles(version):

    # Change version in file Makefile
    change_version_in_file("Makefile", version, "#")

    # Change version in file Makefile.inc
    change_version_in_file("Makefile.inc", version, "#")

    return 2


def change_version_in_fortran_source_files(version):

    rootdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    directories = [os.path.join("src", "adjoint"),
                   os.path.join("src", "base"),
                   os.path.join("src", "common"),
                   os.path.join("src", "sw_mono")]

    count = 0
    for directory in directories:
        files = glob.glob(os.path.join(rootdir, directory, "*.f90"))
        for file in files:
            fname = os.path.relpath(file, rootdir)
            change_version_in_file(fname, version, "!")
            count += 1

    # # Change version in files in directory src/base
    # files = glob.glob(os.path.join(rootdir, "src", "base", "*.f90"))
    # for file in files:
    #     fname = os.path.relpath(file, rootdir)
    #     change_version_in_file(fname, version, "!")
    #     count += 1

    # # Change version in files in directory src/common
    # files = glob.glob(os.path.join(rootdir, "src", "common", "*.f90"))
    # for file in files:
    #     fname = os.path.relpath(file, rootdir)
    #     change_version_in_file(fname, version, "!")
    #     count += 1

    # # Change version in files in directory src/common
    # files = glob.glob(os.path.join(rootdir, "src", "common", "*.f90"))
    # for file in files:
    #     fname = os.path.relpath(file, rootdir)
    #     change_version_in_file(fname, version, "!")
    #     count += 1

    return count


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Modify version in all source files")
    parser.add_argument("version", type=str, help="New version (format:X.X)")
    args = parser.parse_args()

    count = 0
    count += change_version_in_makefiles(args.version)
    count += change_version_in_fortran_source_files(args.version)

    print("MODIFIED FILES: %i" % count)