def gen():
    data = open("inputs/com-dblp/com-dblp.edges")
    output = open("inputs/com-dblp/com-dblp.txt", 'w')
    for line in data.readlines():
        x, y = line[:-1].split("\t")
        if int(x) != int(y):
            output.write(x + " " + y + "\n")
            output.write(y + " " + x + "\n")
    data.close()
    output.close()


if __name__ == "__main__":
    gen()
