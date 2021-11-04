import netCDF4 as nc
import glob
import numpy as np

def main():
    Input_folder = 'MERRA-2/'
    Output_folder = 'MERRA-2/'

    ##### 处理flux nc数据
    data_list_flx = glob.glob(Input_folder + 'MERRA2_400.tavg1_2d_flx*.nc4')
    QLMLs = []
    SPEEDs = []
    TLMLs = []

    for i in range(len(data_list_flx)):
        data = data_list_flx[i]
        nc_data_obj = nc.Dataset(data)
        print(nc_data_obj.variables)
        read_time = np.array(nc_data_obj.variables['time'])  # time, dimension:24
        read_lon = np.array(nc_data_obj.variables['lon'])    # longitude, dimension:576
        read_lat = np.array(nc_data_obj.variables['lat'])    # latitude, dimension:361
        read_QLML = np.array(nc_data_obj.variables['QLML'])        # specific humidity, dimension:(24, 361, 576)
        read_SPEED = np.array(nc_data_obj.variables['SPEED'])      # wind speed, dimension:(24, 361, 576)
        read_TLML = np.array(nc_data_obj.variables['TLML'])      # air temperature, dimension:(24, 576)

        QLMLs.append( np.mean(read_QLML[9: 18, :, :], axis=0) )
        SPEEDs.append( np.mean(read_SPEED[9: 18, :, :], axis=0) )
        TLMLs.append( np.mean(read_TLML[9: 18, :, :], axis=0) )

    QLML = np.mean(QLMLs, axis=0)
    SPEED = np.mean(SPEEDs, axis=0)
    TLML = np.mean(TLMLs, axis=0)

    ##### asm nc数据
    data_list_flx = glob.glob(Input_folder + 'MERRA2_400.tavg3_3d_asm*.nc4')
    RHs = []
    ATs = []

    for i in range(len(data_list_flx)):
        data = data_list_flx[i]
        nc_data_obj = nc.Dataset(data)
        read_RH = np.array(nc_data_obj.variables['RH'])  # relative humidity after moist, dimension:(8, 72, 361, 576)
        read_AT = np.array(nc_data_obj.variables['T'])

        RHs.append(np.mean(read_RH[4: 5, 71, :, :], axis=0))
        ATs.append(np.mean(read_AT[4: 5, 71, :, :], axis=0))

    RHmin = np.min(RHs, axis=0)
    RHmax = np.max(RHs, axis=0)
    ATmin = np.min(ATs, axis=0)-273.15
    ATmax = np.max(ATs, axis=0)-273.15

    ETmin = 0.61*np.exp( 17.27*(ATmin)/(ATmin+237.3) )
    ETmax = 0.61 * np.exp(17.27 * (ATmax) / (ATmax+ 237.3))
    EA = (ETmin*RHmax+ETmax*RHmin)*5

    #写入nc文件
    flx_nc = nc.Dataset(Output_folder + 'flx_rad_0514_0610.nc', 'w', format='NETCDF4')  # 创建nc
    # 确定基础变量的维度信息。相对与坐标系的各个轴(x,y,z)
    flx_nc.createDimension('time', 4)
    flx_nc.createDimension('lat', 361)
    flx_nc.createDimension('lon', 576)

    # 创建变量。参数依次为：‘变量名称’，‘数据类型’，‘基础维度信息’
    flx_nc.createVariable('time', np.int, ('time'))
    flx_nc.createVariable('lat', np.float32, ('lat'))
    flx_nc.createVariable('lon', np.float32, ('lon'))

    # 写入变量time, lat, lon的数据
    time = np.array([5, 14, 6, 10])
    flx_nc.variables['time'][:] = time
    flx_nc.variables['lon'][:] = read_lon
    flx_nc.variables['lat'][:] = read_lat

    # 创建多维变量
    flx_nc.createVariable('QLML', np.float32, ('lat', 'lon'))
    flx_nc.variables['QLML'][:] = QLML

    flx_nc.createVariable('SPEED', np.float32, ('lat', 'lon'))
    flx_nc.variables['SPEED'][:] = SPEED

    flx_nc.createVariable('TLML', np.float32, ('lat', 'lon'))
    flx_nc.variables['TLML'][:] = TLML

    flx_nc.createVariable('VP', np.float32, ('lat', 'lon'))
    flx_nc.variables['VP'][:] = EA

    print(str(len(QLMLs))+' flux files accomplished!')

    ##### 处理radiance nc数据
    data_list_rad = glob.glob(Input_folder + 'MERRA2_400.tavg1_2d_rad*.nc4')
    LWGABs = []
    SWGDNs = []

    for i in range(len(data_list_rad)):
        data = data_list_rad[i]
        nc_data_obj = nc.Dataset(data)
        # print(nc_data_obj.variables)
        read_LWGAB = np.array(nc_data_obj.variables['LWGAB'])  # surface absorbed longwave radiation:(24, 361, 576)
        read_SWGDN = np.array(nc_data_obj.variables['SWGDN'])  # surface incoming shortwave flux, dimension:(24, 361, 576)

        LWGABs.append( np.mean(read_LWGAB[9: 18, :, :], axis=0) )
        SWGDNs.append( np.mean(read_SWGDN[9: 18, :, :], axis=0) )

    LWGAB = np.mean(LWGABs, axis=0)
    SWGDN = np.mean(SWGDNs, axis=0)

    # 创建多维变量
    flx_nc.createVariable('LWGAB', np.float32, ('lat', 'lon'))
    flx_nc.variables['LWGAB'][:] = LWGAB

    flx_nc.createVariable('SWGDN', np.float32, ('lat', 'lon'))
    flx_nc.variables['SWGDN'][:] = SWGDN

    # 关闭文件
    flx_nc.close()
    print(str(len(SWGDNs))+ ' radiance files accomplished!')
    print('Success!')


main()


########################################################################

