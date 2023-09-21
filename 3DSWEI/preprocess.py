import os
import cv2 as cv


def remove9999(start=1, end=9999):
    fileRt = path + 'fileRt.txt'
    fileRtBack = path + 'fileRtBack.txt'

    # imgPath = path + "USIMG\\"
    # imgName = "USImage"
    imgFormat = "bmp"

    imgPath = path + "test\\"
    # imgPath1 = path + "test - 副本 (2)\\"
    imgName = ""
    imgFormat = "mat"

    recnt = 0
    with open(fileRt, mode='r') as iFile:
        with open(fileRtBack, mode='w') as oFile:
            for order in range(start, end + 1):
                line = iFile.readline()
                if '-9999' in line:
                    iFile.readline()
                    recnt += 1
                    os.remove(f'{imgPath}{imgName}{order}.{imgFormat}')
                    # os.remove(f'{imgPath1}{imgName}{order}.{imgFormat}')
                    continue
                oFile.write(line)
                line = iFile.readline()
                oFile.write(line)

                os.rename(f'{imgPath}{imgName}{order}.{imgFormat}',
                          f'{imgPath}{imgName}{order - recnt}.{imgFormat}')


def cropPicture(top, bot, left, right, start=1, end=999):
    imgPath = path + "USIMG\\"
    imgBackPath = path + "USIMG\\"
    imgName = "USImage"
    imgFormat = "bmp"
    for i in range(start, end + 1):
        img = cv.imread(f'{imgPath}{imgName}{i}.{imgFormat}')
        img = img[top:bot + 1, left:right + 1]
        cv.imwrite(f'{imgBackPath}{imgName}{i}.{imgFormat}', img)


def transcolor(start, end):
    readpath = path + 'USIMG\\'
    savepath = path + 'USIMG\\'
    imgName = "USImage"
    imgFormat = "bmp"
    for i in range(start, end + 1):
        img = cv.imread(f'{readpath}{imgName}{i}.{imgFormat}', cv.IMREAD_GRAYSCALE)
        cv.normalize(img, img, 0, 255, cv.NORM_MINMAX)
        cv.imwrite(f'{savepath}{imgName}{i}.{imgFormat}', img)


def renameMat(start, end):
    imgPath = path + "test\\"
    imgName = ""
    imgFormat = "mat"
    recnt = 2
    for order in range(start, end + 1):
        os.rename(f'{imgPath}{imgName}{order}.{imgFormat}',
                  f'{imgPath}{imgName}{order - recnt}.{imgFormat}')
    print(f'重命名.mat文件，从{start}.mat到{end}.mat，修改为{start - recnt}.mat到{end - recnt}.mat')


if __name__ == '__main__':
    # path = 'Recon' + '1122\\'
    path = 'E:\\nju307_wt\\SWI_m\\1208\\1\\'
    # filenum = len(os.listdir(path + 'USIMG'))
    filenum = len(os.listdir(path + 'test'))
    # renameMat(start=3, end=filenum + 2)
    remove9999(start=1, end=filenum)
    # cropPicture(top=21, bot=303, left=53, right=609, start=1, end=filenum)
    # transcolor(start=1, end=filenum)
    # if not os.path.exists(path + "ReconsDataP"):
    #     os.makedirs(path + "ReconsDataP")
