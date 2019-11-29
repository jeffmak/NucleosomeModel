
with open("analyzed_data_complete_bw20.txt") as in_file:
    content = in_file.read()
    contents = content.split('# S')
    for content in contents:
        with open('S'+content.split('\n')[0], 'w+') as w_file:
            w_file.write('S'+content)