import csv

def filter_nodes(file_path,nodes_path):
    with open(nodes_path, 'r') as f:
        nodes = set(line.strip() for line in f)

    with open(file_path, 'r') as f, open(f'reduced_{file_path}', 'w') as out:
        reader = csv.reader(f, delimiter=' ')
        writer = csv.writer(out, delimiter=' ')
        for row in reader:
            if row[0] in nodes and row[1] in nodes:
                writer.writerow(row)