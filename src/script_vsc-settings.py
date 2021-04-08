# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 14:30:33 2021

@author: linigodelacruz
"""
#!/usr/bin/env python3
from os import path
import tempfile
import subprocess
import datetime as dt
import click


@click.group()
def cli():
    pass


@cli.command(name='export')
@click.option('--settings-dir', default=path.expanduser('~/.config/Code/User'))
@click.option('--out-dir', default='./')
def export_settings(settings_dir, out_dir):
    out_dir = path.abspath(out_dir)
    dump_path = tempfile.mkdtemp()
    click.echo('Exporting vsc settings:')
    click.echo('created a temporary dump dir {}'.format(dump_path))
    ext_file = path.join(dump_path, 'extensions.txt')
    click.echo('generating extensions list')
    subprocess.run('code --list-extensions > {}'.format(ext_file), shell=True)

    for f in ['settings.json', 'keybindings.json', 'projects.json', 'snippets']:
        ff = path.join(settings_dir, f)
        if path.exists(ff):
            click.echo('copying ' + ff)
            subprocess.run('cp -r {} {}'.format(ff, dump_path), shell=True)

    now = dt.datetime.now()
    ball_path = path.join(out_dir, 'vsc-settings-{}.zip'.format(now.strftime("%Y-%m-%d-%H%M%S")))
    subprocess.run(
        'cd {} && zip -r {} .'.format(dump_path, ball_path),
        shell=True
    )
    click.echo('VSC settings exported into ' + ball_path)


@cli.command(name='import')
def import_settings():
    click.echo('Importing settings.')


if __name__ == "__main__":
    cli()