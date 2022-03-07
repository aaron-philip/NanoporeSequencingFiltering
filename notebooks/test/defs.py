from multiprocessing import Queue, Process
def add_helper(queue, arg1, arg2): # the func called in child processes
    ret = arg1 + arg2
    queue.put(ret)
 
def multi_add(): # spawns child processes
    q = Queue()
    processes = []
    rets = []
    for _ in range(0, 100):
        p = Process(target=add_helper, args=(q, 1, 2))
        processes.append(p)
        p.start()
    for p in processes:
        ret = q.get() # will block
        rets.append(ret)
    for p in processes:
        p.join()
    return rets

def returnFirst():
    output = multi_add()
    return(output[0])

if __name__=="__main__":
    print(returnFirst())
