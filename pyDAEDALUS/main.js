const { app, BrowserWindow } = require('electron');
const {ipcMain, dialog} = require('electron')
const path = require ('path');
const fs = require('fs');
const os = require('os');
// var zerorpc = require("zerorpc");

// Keep a global reference of the window object, if you don't, the window will
// be closed automatically when the JavaScript object is garbage collected.
let win

function createWindow () {
  // call python DAEDALUS Server
  // var subpy = require('child_process').spawn('./dist/DAEDserve');
  // Create the window.
  win = new BrowserWindow({
    width: 500,
    height: 440,
    webPreferences: {
      nodeIntegration: true
    }
  })

  // and load the index.html of the app.
  win.loadFile('index.html')

  // Open the DevTools.
  // win.webContents.openDevTools()

  // Emitted when the window is closed.
  win.on('closed', () => {
    // Dereference the window object, usually you would store windows
    // in an array if your app supports multi windows, this is the time
    // when you should delete the corresponding element.
//    subpy.kill('SIGINT');
    win = null
  })
}
var plyfile='tet.ply'
var seqfile='M13.txt'
// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.on('ready', createWindow)

// Quit when all windows are closed.
app.on('window-all-closed', () => {
  // On macOS it is common for applications and their menu bar
  // to stay active until the user quits explicitly with Cmd + Q
  if (process.platform !== 'darwin') {
    app.quit()
  }
})

app.on('activate', () => {
  // On macOS it's common to re-create a window in the app when the
  // dock icon is clicked and there are no other windows open.
  if (win === null) {
    createWindow()
  }
})

// In this file you can include the rest of your app's specific main process
// code. You can also put them in separate files and require them here.

ipcMain.on('open-plyfile-dialog', (event) => {
  let plyfiles = dialog.showOpenDialogSync({ properties: ['openFile'] })
  plyfile = plyfiles[0]
  win.webContents.send('plyTableUpdate',plyfile);
})
ipcMain.on('open-seqfile-dialog', (event) => {
  let seqfiles = dialog.showOpenDialogSync({ properties: ['openFile'] })
  seqfile = seqfiles[0]
  win.webContents.send('seqTableUpdate',seqfile);
})
ipcMain.on('changeDX', (event, args) => {
  dxType = args
  win.webContents.send('dxUpdate',dxType);
})

exports.handleForm = function handleForm(targetWindow, projName, helicalForm, helicalTurns) {
//  var client = new zerorpc.Client();
//  client.connect("tcp://127.0.0.1:4242");
  
  win.webContents.send('finishedUpdate',"Processing...", projName, helicalForm, helicalTurns, plyfile, seqfile);
//  client.invoke("calc", projName, helicalForm, helicalTurns, plyfile, seqfile, function(error, results) {
//    console.log(error);
//    let tmpSt=results.toString()
//    console.log(tmpSt);
//    if (tmpSt==="Finished!") {
//      win.webContents.send('finishedUpdate',"Success!");
//    }
//  });
};
