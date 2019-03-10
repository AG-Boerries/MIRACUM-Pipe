MIRACUM workflows for Galaxy
============================

Users
-----

These are the Galaxy Workflow files that define the MIRACUM pipeline used in
Galaxy. To use these workflows, import the files into any Galaxy instance that
has all required tools installed (currrently only https://usegalaxy.eu), then
run the *Miracum - main*, *Miracum - variant annotation* and *Miracum - report
variants* workflows, in this order, to go from raw sequenced exome reads to
annotated tabular variant reports.


Developers
----------

Galaxy Workflow files are JSON-formatted. You can edit them manually, or import
them into Galaxy and use Galaxy's graphical workflow editor, then download the
edited version again. Either way, you should pretty-format the edited file
before committing them to obtain nice and readable diffs.
Use Python's built-in ``json.tool`` command for this task::

  python3 -m json.tool your_edited_workflow.ga > workflow_to_commit.ga
  
Python3.5 or later is required because earlier versions won't preserve element
order.

When you edit a subworkflow, please note that the changes you make do not get
reflected in the containing workflow (which carries its own embedded copy of
the subworkflow). To make your changes visible in the containing workflow, you
need to remove the embedded subworkflow from it and add the edited version back
in.

