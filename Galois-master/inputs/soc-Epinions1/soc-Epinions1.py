def gen_epinions():
    data = open("inputs/soc-Epinions1/soc-Epinions1.edges")
    output = open("inputs/soc-Epinions1/soc-Epinions1.txt", 'w')
    for line in data.readlines():
        x, y = line[:-1].split("\t")
        if int(x) != int(y):
            output.write(x + " " + y + "\n")
            output.write(y + " " + x + "\n")
    data.close()
    output.close()


if __name__ == "__main__":
    gen_epinions()
