import matplotlib.pyplot as plt

def plot(x_list, y_list):
    plt.plot(x_list, y_list)
    plt.show()

job = Job.query.filter(Job.id==2).first()
result = json.loads(job.result)
plot(result['t_list'], result['vis_list'])