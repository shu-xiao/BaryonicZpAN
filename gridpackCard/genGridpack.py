import os
import shutil

def mkNdir(dirName):
    if not os.path.isdir(dirName): os.mkdir(dirName)
def getWidth(mzp, gq=0.25):
    width = 0
    f_banner = open("cards/ZpBaryonic_bb_bannerDir/ZpBaryonic_bb_MZp" + str(mzp) + "_gq" + str(format(100*gq,'.0f')) + "_MDM100_banner.txt")
    for line in f_banner:
        if line.find("900000") > 0 and line.find("DECAY") >=0:
            num0 = line.split("9000001")[1].split('e')[0]
            num1 = line.split("9000001")[1].split('e')[1]
            width = float(num0)*(10**float(num1))
            break
    print(width)
    return width

def main():
    gq=0.25
    #for mzp in range(1500,1600,100):
    for mzp in range(600,2100,100):
        if (mzp!=1000 and mzp!=1100 and mzp!=600): continue
        dirName = 'ZpBaryonic_bb_MZp'+str(mzp)+'_gq'+str(format(100*gq,'.0f'))+'_MDM100'
        mkNdir('cards/'+dirName)
        print('create '+dirName)
        shutil.copyfile('cards/ZpBaryonic_bb_template/tem_run_card.dat','cards/'+dirName+'/'+dirName+'_run_card.dat')
        shutil.copyfile('cards/ZpBaryonic_bb_template/tem_extramodels.dat','cards/'+dirName+'/'+dirName+'_extramodels.dat')

        f_proc0 = open('cards/ZpBaryonic_bb_template/tem_proc_card.dat','r')
        f_proc1 = open('cards/'+dirName+'/'+dirName+'_proc_card.dat','w')
        for line in f_proc0:
            f_proc1.write(line.replace('ZpBaryonic_bb_MZp1000_gq25_MDM100',dirName))
        f_proc0.close()
        f_proc1.close()

        f_cust0 = open('cards/ZpBaryonic_bb_template/tem_customizecards.dat','r')
        f_cust1 = open('cards/'+dirName+'/'+dirName+'_customizecards.dat','w')
        for line in f_cust0:
            if line.find('MZP') > 0: f_cust1.write(line.replace('MZP',str(mzp)))
            elif line.find('GZ') > 0: f_cust1.write(line.replace('GZ',str(mzp)))
            #elif line.find('Width') > 0: f_cust1.write(line.replace('Width',str(84.5384)))
            elif line.find('Width') > 0: f_cust1.write(line.replace('Width',str(format(getWidth(mzp), '.4f'))))
            else: f_cust1.write(line)
        f_cust0.close()
        f_cust1.close()
        command = './gridpack_generation.sh ' + dirName + ' cards/' + dirName + ' 1nd'
        print(command)
        os.system(command)

if __name__ == "__main__":
    main()
    #os.system('cp *.tarball.tar.xz ~/public/ZpBaryonic_gridpack_new')
