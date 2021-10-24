def FolderCreater(path):
    '''
    Create a new folder    
    path: for example (path = r'D:\Example\Folder_1')
    '''
    import os
    if not os.path.exists(path):
        os.makedirs(path)
