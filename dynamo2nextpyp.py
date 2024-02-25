#!/usr/bin/env python3

import sys
import os
import numpy as np
from eulerangles import ConversionMeta, convert_eulers

def parse_dynamo_tbl(file_path):
    data = []
    tomogram_particle_counts = {}

    with open(file_path, 'r') as file:
        for line in file:
            values = line.strip().split()

            if len(values) > 19:
                tomogram_number = int(values[19])

                if tomogram_number not in tomogram_particle_counts:
                    tomogram_particle_counts[tomogram_number] = 0
                else:
                    tomogram_particle_counts[tomogram_number] += 1
                
                row_data = {
                    'tomo': tomogram_number,
                    'tdrot': float(values[6]),
                    'tilt': float(values[7]),
                    'narot': float(values[8]),
                    'particle_index': tomogram_particle_counts[tomogram_number]
                }

                data.append(row_data)
            else:
                print(f"Warning: Row with insufficient columns found: {line}")

    return data

def parse_columns20_index(file_path):
    index_map = {}

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            tomogram_index = int(parts[0])
            full_filename = ' '.join(parts[1:])
            basename = os.path.splitext(os.path.basename(full_filename))[0]
            index_map[tomogram_index] = basename

    return index_map

def parse_volumes_list(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    header = lines[0].strip().split()
    data = []

    for line in lines[1:]:
        values = line.strip().split()
        row_data = {'original_line': line.strip()}
        for i, h in enumerate(header):
            if h == 'filename':
                filename = values[i]
                basename = os.path.splitext(os.path.basename(filename))[0]
                parts = basename.split('_spk')
                tomogram_name = parts[0] if parts else 'unknown'
                particle_index_str = parts[1] if len(parts) > 1 else '0'
                row_data['tomogram_name'] = tomogram_name
                row_data['particle_index'] = int(particle_index_str)
            else:
                try:
                    row_data[h] = float(values[i]) if '.' in values[i] else int(values[i])
                except ValueError:
                    row_data[h] = values[i]
        data.append(row_data)

    return data, header

def write_updated_volumes_list(data, output_file, header):
    with open(output_file, 'w') as file:
        file.write('\t'.join(header) + '\n')

        for row in data:
            original_values = row['original_line'].split()
            normal_indices = [header.index('normalX'), header.index('normalY'), header.index('normalZ')]
            for i, normal_index in enumerate(normal_indices):
                original_values[normal_index] = str(row[f'normalX' if i == 0 else f'normalY' if i == 1 else 'normalZ'])
            file.write('\t'.join(original_values) + '\n')

def update_volumes_list_with_euler_angles(dynamo_file, columns20_file, volumes_list_file, output_file):
    dynamo_data = parse_dynamo_tbl(dynamo_file)
    columns20_map = parse_columns20_index(columns20_file)
    volumes_list_data, header = parse_volumes_list(volumes_list_file)

    source_metadata = ConversionMeta(name='input', axes='zxz', intrinsic=False, right_handed_rotation=True, active=False)
    target_metadata = ConversionMeta(name='output', axes='zxz', intrinsic=False, right_handed_rotation=False, active=False)

    # Define and populate euler_map
    euler_map = {}
    for row in dynamo_data:
        tomogram_name = columns20_map.get(row['tomo'])
        if tomogram_name:
            key = (tomogram_name, row['particle_index'])
            euler_angles = ([row['tdrot'], row['tilt'], row['narot']])
            euler_map[key] = euler_angles

    for row in volumes_list_data:
        tomogram_name_from_list = row['tomogram_name']
        
        key = (row['tomogram_name'], row['particle_index'])
        if key in euler_map:
            input_eulers = np.array(euler_map[key]).reshape((1, 3))
         
            output_eulers = convert_eulers(input_eulers, source_meta=source_metadata, target_meta=target_metadata)
            alpha, beta, gamma = output_eulers.flatten()
            normalZ = -alpha
            normalX = -beta
            normalY = -gamma
            
            row['normalX'], row['normalY'], row['normalZ'] = normalX, normalY, normalZ

    write_updated_volumes_list(volumes_list_data, output_file, header)

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python script.py dynamo_file columns20_file volumes_list_file output_file")
        sys.exit(1)

    dynamo_file, columns20_file, volumes_list_file, output_file = sys.argv[1:]
    update_volumes_list_with_euler_angles(dynamo_file, columns20_file, volumes_list_file, output_file)

