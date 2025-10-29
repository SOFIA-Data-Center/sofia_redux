#!/bin/sh
pyuic6 main.ui -o ../ui/ui_main.py
pyuic6 pipe_step.ui -o ../ui/ui_pipe_step.py
pyuic6 textview.ui -o ../ui/ui_textview.py
pyuic6 qad_settings.ui -o ../ui/ui_qad_settings.py
pyuic6 remove_files.ui -o ../ui/ui_remove_files.py
pyuic6 edit_param.ui -o ../ui/ui_edit_param.py
pyuic6 progress.ui -o ../ui/ui_progress.py
