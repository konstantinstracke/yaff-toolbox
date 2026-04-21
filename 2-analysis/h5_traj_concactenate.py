import h5py
import numpy as np
import sys

def map_integers_to_strings(ids, strings):
    return [strings[i].decode() if isinstance(strings[i], bytes) else strings[i] for i in ids]

def read_system_group(f):
    system_data = {}
    grp = f['system']
    for key in grp:
        system_data[key] = np.array(grp[key])
    return system_data

def read_trajectory_group(f):
    traj_data = {}
    grp = f['trajectory']
    for key in grp:
        traj_data[key] = np.array(grp[key])
    return traj_data

def concat_trajectories(h5_files, output_file):
    all_traj = {}
    traj_shapes = {}
    system_data = None

    for i, hfile in enumerate(h5_files):
        print(f"Reading: {hfile}")
        with h5py.File(hfile, 'r') as f:
            if i == 0:
                system_data = read_system_group(f)

            traj = read_trajectory_group(f)

            for key, arr in traj.items():
                traj_shape = arr.shape[1:]  # Shape excluding frame count
                if key in traj_shapes:
                    if traj_shapes[key] != traj_shape:
                        raise ValueError(
                            f"Shape mismatch for '{key}': expected shape (*, {traj_shapes[key]}), "
                            f"but got (*, {traj_shape}) in file {hfile}"
                        )
                else:
                    traj_shapes[key] = traj_shape

                if key not in all_traj:
                    all_traj[key] = [arr]
                else:
                    all_traj[key].append(arr)

    for key in all_traj:
        all_traj[key] = np.vstack(all_traj[key])

    print(f"Writing to: {output_file}")
    with h5py.File(output_file, 'w') as f_out:
        sys_grp = f_out.create_group('system')
        for key, value in system_data.items():
            dtype = 'S22' if value.dtype.kind in {'U', 'S'} else None
            sys_grp.create_dataset(key, data=value.astype(dtype) if dtype else value)

        traj_grp = f_out.create_group('trajectory')
        for key, value in all_traj.items():
            traj_grp.create_dataset(key, data=value)

    print("Done.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python h5_concatenate.py output_file.h5 input1.h5 input2.h5 ...")
        print("Example: python h5_concactenate.py traj_all.h5 traj_{0..57}.h5")
        print("Same data shapes needed !")
        sys.exit(1)

    output_file = sys.argv[1]
    input_files = sys.argv[2:]

    concat_trajectories(input_files, output_file)
