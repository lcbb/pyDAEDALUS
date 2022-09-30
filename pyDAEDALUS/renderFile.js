const { remote, ipcRenderer } = require('electron');
const { handleForm} = remote.require('./main');
const currentWindow = remote.getCurrentWindow();

const submitFormButton = document.querySelector("#ipcForm2");
const responseParagraph = document.getElementById('response');
const dximgParagraph = document.getElementById('DXimg');
const plyTableCell = document.getElementById('plyTable');
const seqTableCell = document.getElementById('seqTable');
const helF = document.getElementById("helicalForm");
const selectPlyBtn = document.getElementById('select-ply');
const selectSeqBtn = document.getElementById('select-seq');

submitFormButton.addEventListener("submit", function(event){
        event.preventDefault();   // stop the form from submitting
        let projName = document.getElementById("projName").value;
        let helicalForm = document.getElementById("helicalForm").value;
        let helicalTurns = document.getElementById("helicalTurns").value;
        handleForm(currentWindow, projName, helicalForm, helicalTurns)
});

ipcRenderer.on('form-received', function(event, args){
    responseParagraph.innerHTML = args
/*
        you could choose to submit the form here after the main process completes
        and use this as a processing step
*/
});

ipcRenderer.on('plyTableUpdate', function(event, args){
  plyTableCell.innerHTML = args
});

ipcRenderer.on('seqTableUpdate', function(event, args){
  seqTableCell.innerHTML = args
});

ipcRenderer.on('finishedUpdate', function(event, args, projNameA, helicalFormA, helicalTurnsA, plyfileA, seqfileA){
  responseParagraph.innerHTML = args
  $.xmlrpc({
    url: 'http://localhost:4242',
    methodName: 'calc',
    params: [projNameA, helicalFormA, helicalTurnsA, plyfileA, seqfileA],
    success: function(response, status, jqXHR) { responseParagraph.innerHTML = "Success!" },
    error: function(jqXHR, status, error) { responseParagraph.innerHTML = "Fail!" }
  });
});

ipcRenderer.on('dxUpdate', function(event, args){
  if (args==="Bform") {
    dximgParagraph.innerHTML = '<p id="DXimg"><img src="./images/DX_Bform.png" height="80"></p>'
  }
  if (args==="Aform") {
    dximgParagraph.innerHTML = '<p id="DXimg"><img src="./images/DX_Aform.png" height="80"></p>'
  }
  if (args==="Hybrid") {
    dximgParagraph.innerHTML = '<p id="DXimg"><img src="./images/DX_Hform.png" height="80"></p>'
  }
  if (args==="Twisted") {
    dximgParagraph.innerHTML = '<p id="DXimg"><img src="./images/DX_Altform.png" height="80"></p>'
  }
});

selectPlyBtn.addEventListener('click', (event) => {
  ipcRenderer.send('open-plyfile-dialog')
})

selectSeqBtn.addEventListener('click', (event) => {
  ipcRenderer.send('open-seqfile-dialog')
})

helF.addEventListener('change', (event) => {
  let dxChoice = document.getElementById("helicalForm").value;
  console.log(dxChoice)
  ipcRenderer.send('changeDX', dxChoice)
});
