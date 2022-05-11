def gen_citeseer():
    data = open("inputs/citeseer/citeseer.edges")
    output = open("inputs/citeseer/citeseer.txt", 'w')
    for line in data.readlines():
        x, y, l = line[:-1].split(",")
        if int(x) != int(y):
            output.write(x + " " + y + "\n")
            output.write(y + " " + x + "\n")
    data.close()
    output.close()


if __name__ == "__main__":
    gen_citeseer()
