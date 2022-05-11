def gen():
    data = open("inputs/com-youtube/com-youtube.edges")
    output = open("inputs/com-youtube/com-youtube.txt", 'w')
    for line in data.readlines():
        x, y = line[:-1].split("\t")
        x = int(x)
        y = int(y)
        if x != y:
            output.write(str(x-1) + " " + str(y-1) + "\n")
            output.write(str(y-1) + " " + str(x-1) + "\n")
    data.close()
    output.close()


if __name__ == "__main__":
    gen()
