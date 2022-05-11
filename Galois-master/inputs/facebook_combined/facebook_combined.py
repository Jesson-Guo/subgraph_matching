def gen_facebook_combined():
    data = open("inputs/facebook_combined/facebook_combined.edges")
    output = open("inputs/facebook_combined/facebook_combined.txt", 'w')
    for line in data.readlines():
        x, y = line[:-1].split(" ")
        if int(x) != int(y):
            output.write(x + " " + y + "\n")
            output.write(y + " " + x + "\n")
    data.close()
    output.close()


if __name__ == "__main__":
    gen_facebook_combined()
