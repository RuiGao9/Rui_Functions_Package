def FolderCreater(path):
    '''
    path: for example (path = r'D:\Example\Folder_1')
            Folder_1 is the new folder
    '''
    import os
    if not os.path.exists(path):
        os.makedirs(path)