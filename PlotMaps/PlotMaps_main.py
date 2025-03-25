import matplotlib.pyplot as plt  
from mpl_toolkits.basemap import Basemap  
def draw_map():  
    map = Basemap(projection='merc', llcrnrlat=-60, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')  
    map.drawcoastlines()  
    map.drawcountries()  
    plt.show()  
if __name__ == "__main__":  
    draw_map()