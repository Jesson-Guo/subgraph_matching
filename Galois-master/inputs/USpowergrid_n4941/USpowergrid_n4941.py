def gen():
    data = open("inputs/USpowergrid_n4941/USpowergrid_n4941.edges")
    output = open("inputs/USpowergrid_n4941/USpowergrid_n4941.txt", 'w')
    for line in data.readlines():
        x, y, _ = line[:-1].split(" ")
        if int(x) != int(y):
            output.write(x + " " + y + "\n")
            output.write(y + " " + x + "\n")
    data.close()
    output.close()


if __name__ == "__main__":
    gen()
