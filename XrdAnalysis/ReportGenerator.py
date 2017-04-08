import logging
import os
import re
import subprocess
import sys

import jinja2

import XrdAnalysis
from XrdAnalysis import Reader


class Generator(object):
    def __init__(self):
        self.context_dict = None

    @staticmethod
    def render(template_path, context):
        path, filename = os.path.split(template_path)
        return jinja2.Environment(
            loader=jinja2.FileSystemLoader(
                path or './')).get_template(filename).render(context)

    def print_to_tex(self, file_list):
        class_path = os.path.dirname(XrdAnalysis.__file__)
        xrd_block = [Reader.reader(i).read_data().print_result_tex()
                     for i in file_list]
        body_dict = {'xrd_block': xrd_block}

        tmp_file = os.path.join(class_path, 'templates', 'tmp.tex')
        logging.info("Use template file {0}".format(os.path.abspath(tmp_file)))
        tex_str = self.render(tmp_file, body_dict)
        tex_str = re.sub(r'\{ ', '{', tex_str)
        tex_str = re.sub(r' \}', '}', tex_str)

        logging.debug(tex_str)

        return tex_str

    def print_to_pdf(self, file_list, **param):
        tex_str = self.print_to_tex(file_list)

        if 'file_name' in param:
            tex_file_name = param['file_name'].replace('.pdf', '.tex')
        else:
            tex_file_name = 'temp.tex'
        with open(tex_file_name, 'wb') as th:
            th.write(tex_str.encode('utf-8'))
        tex_file = os.path.abspath(tex_file_name)
        process = subprocess.Popen(
            'pdflatex --interaction=nonstopmode {0}'.format(tex_file),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        out, err = process.communicate()
        # if out:
        #     print("standard output of subprocess:")
        #     print(out)
        # if err:
        #     print("standard error of subprocess:")
        #     print(err)


if __name__ == '__main__':
    logging.basicConfig(
        # filename=os.path.join(
        #     os.path.dirname(sys.argv[0]), 'log', __name__ + '.log'),
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
    )
    s = Generator()
    s.print_to_pdf([os.path.join(
        os.path.dirname(sys.argv[0]),
        'TestData', 'pf.raw'
    )])
