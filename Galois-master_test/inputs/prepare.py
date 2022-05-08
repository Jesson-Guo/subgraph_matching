def gen_complete(size):
    complete_edges = open("inputs/wiki-vote/complete.txt", "w")
    for i in range(size+1):
        for j in range(size+1):
            if i==j:
                continue
            else:
                complete_edges.write(str(i) + " " + str(j) + "\n")
    complete_edges.close()


def gen_dataset(path):
    data = open(path)
    output = open("test/test2.txt", 'w')
    for line in data.readlines():
        x, y = line[:-1].split(" ")
        if int(x) != int(y):
            output.write(x + " " + y + "\n")
            output.write(y + " " + x + "\n")
    data.close()
    output.close()


if __name__ == "__main__":
    gen_dataset("/home/jesson/Documents/code/GraphPi-master/dataset/test_input")
