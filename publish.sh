echo -e "
A script that builds whl
package and uploads it.

Version: v1.0.0
Authors:
	name: Platon Bykadorov
	email: platon.work@gmail.com
	years: 2024\n"

cd $(dirname $(realpath $0))
python3 -m pip install --upgrade build; echo
python3 -m pip install --upgrade twine; echo
if [[ -d dist ]]; then
	rm -vr dist; echo
fi
python3 -m build; echo
python3 -m twine upload dist/*; echo
rm -vr dist; echo