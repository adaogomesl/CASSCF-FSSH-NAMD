import sys

def getindex(index):
    ## This function read single, range, separete range index and convert them to a list
    index_list=[]
    for i in index:
        if '-' in i:
            a,b=i.split('-')
            a,b=int(a),int(b)
            index_list+=range(a,b+1)
        else:
            index_list.append(int(i))    

    index_list=sorted(list(set(index_list))) # remove duplicates and sort from low to high
    return index_list

def main():
    if len(sys.argv) < 1:
        print(info)
        exit()

    with open('range-s0.txt','r') as para:
        follow=para.read().split()
    follow=getindex(follow)
    print(follow)
    print(*follow, sep=' ')
if __name__ == '__main__':
    main()

